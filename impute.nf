params.pedFile = 'TODO'
params.mapFile = 'TODO'
params.genticMapDir = 'TODO'
params.chromosomeSizesFile = 'b37.chrom.sizes'
params.referenceHapsFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.haplotypes.gz"
params.referenceLegendFilePattern = "ALL.chr%s.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.nomono.legend.gz"
params.referenceGeneticMapPattern = "genetic_map_chr%s_combined_b37.txt"
params.referenceSample = "ALL.integrated_phase1_v3.20101123.snps_indels_svs.genotypes.sample"
params.publishDirPath = "$HOME/imputation-results" 

mapFileChan = Channel.fromPath(params.mapFile)
pedFileChan = Channel.fromPath(params.pedFile)
chromosomesList = 1..22

db_path = file(params.genticMapDir)

process plink {

  container 'insilicodb/docker-impute2'

  maxForks 6

  input:
  file mapFile from mapFileChan
  file pedFile from pedFileChan
  each chromosome from chromosomesList 

  output:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into plinkOutChan

  """
  plink-1.9 --noweb --file ${pedFile.baseName} --chr $chromosome --make-bed --out chr${chromosome}
  plink-1.9 --bfile chr${chromosome} --list-duplicate-vars ids-only suppress-first
  [[ -e "plink.dupvar" ]] && plink-1.9 --bfile chr${chromosome} --exclude plink.dupvar --make-bed --out chr${chromosome}
  """

}

//plinkOutChan.subscribe{  println "${it.name}"}


process shapeitCheck {
  validExitStatus 0,1,2
  errorStrategy 'ignore'

  maxForks 6

  container 'insilicodb/docker-impute2'

  input:
  set val(chromosome), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") from plinkOutChan
  file db_path

  output:
  set val(chromosome), file("chr${chromosome}.alignments.log"), file("chr${chromosome}.alignments.snp.strand.exclude"), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") into shapitCheckChan

  script:
  hapFile = file( db_path.name + "/" + sprintf(params.referenceHapsFilePattern, chromosome) )
  legendFile = file( db_path.name + "/" + sprintf(params.referenceLegendFilePattern, chromosome) )
  sampleFile = file( db_path.name + "/" + params.referenceSample )

  """
  shapeit -check --input-bed chr${chromosome}.bed chr${chromosome}.bim chr${chromosome}.fam --input-ref $hapFile $legendFile $sampleFile --output-log chr${chromosome}.alignments
  """

}


//shapitCheckChan.subscribe{  println "${it.name}"}

process shapeit {

  maxForks 4

  container 'insilicodb/docker-impute2'

  input:
  set val(chromosome), file("chr${chromosome}.alignments.log"), file("chr${chromosome}.alignments.snp.strand.exclude"), file("chr${chromosome}.bed"), file("chr${chromosome}.fam"), file("chr${chromosome}.bim") from shapitCheckChan 
  file db_path

  output:
  set val(chromosome), file("chr${chromosome}.phased.haps"), file("chr${chromosome}.phased.sample") into shapeitChan

  script:
  hapFile = file( db_path.name + "/" + sprintf(params.referenceHapsFilePattern, chromosome) )
  legendFile = file( db_path.name + "/" + sprintf(params.referenceLegendFilePattern, chromosome) )
  sampleFile = file( db_path.name + "/" + params.referenceSample )
  excludeFile = "chr${chromosome}.alignments.snp.strand.exclude"
  mapFile = file( db_path.name + "/" + sprintf(params.referenceGeneticMapPattern, chromosome) )

  """
  shapeit --input-bed chr${chromosome}.bed chr${chromosome}.bim chr${chromosome}.fam --input-ref $hapFile $legendFile $sampleFile --exclude-snp $excludeFile --input-map $mapFile -O chr${chromosome}.phased 
  """

}

imputeChromChunckChannel = shapeitChan.flatMap { chromosome, hapsFile, sampleFile ->
   def results = []
   
   def chunks = getChromosomeChunkPairs(getChromosomeSize(file(params.chromosomeSizesFile), chromosome))
   
   chunks.each { chunkStart, chunkEnd -> 
     results.push( [ chromosome, hapsFile, sampleFile, chunkStart, chunkEnd] )
   }
   
   return results 
}

process impute2 {

  errorStrategy 'ignore'

  container 'insilicodb/docker-impute2'

  maxForks 2
  
  input:
  set val(chromosome), file("chr${chromosome}.phased.haps"), file("chr${chromosome}.phased.sample"), val(chunkStart), val(chunkEnd) from imputeChromChunckChannel
  file db_path

  output:
  set val(chromosome), file("chr${chromosome}-${chunkStart}-${chunkEnd}.imputed") into impute2Chan

  script:
  hapFile = file( db_path.name + "/" + sprintf(params.referenceHapsFilePattern, chromosome) )
  legendFile = file( db_path.name + "/" + sprintf(params.referenceLegendFilePattern, chromosome) )
  sampleFile = file( db_path.name + "/" + params.referenceSample )
  mapFile = file( db_path.name + "/" + sprintf(params.referenceGeneticMapPattern, chromosome) )

  """
  impute2 -use_prephased_g -known_haps_g chr${chromosome}.phased.haps -h $hapFile -l $legendFile -m $mapFile -int $chunkStart $chunkEnd -Ne 20000 -o chr${chromosome}-${chunkStart}-${chunkEnd}.imputed
  
  #sometimes there are no SNP's in a region
  if [ ! -f "chr${chromosome}-${chunkStart}-${chunkEnd}.imputed" ]; then
    touch "chr${chromosome}-${chunkStart}-${chunkEnd}.imputed";
  fi
  
  """
}

impute2List = impute2Chan.toSortedList() //gives a dataFlow instance, nee to get the val property of it
impute2List = impute2List.val

impute2Map = [:]

impute2List.each { chrom, file ->

  if ( !impute2Map.containsKey(chrom) ) {
   impute2Map.put(chrom, [])
  }
  impute2Map.get(chrom).add(file)

}

impute2MapChannel = Channel.create()

impute2Map.each { chrom, fileList ->
  impute2MapChannel.bind([chrom, fileList])
}

impute2MapChannel.close()

process impute2Concat {
  
  publishDir params.publishDirPath
 
  input:
  set val(chromosome), file(imputedFiles) from impute2MapChannel

  output:
  set val(chromosome), file("chr${chromosome}.imputed") into impute2ConcatChan

  """
  cat $imputedFiles > chr${chromosome}.imputed
  """
}


// ----==== utility methods ====----


def getChromosomeSize( chromosomeSizesFile, chromosome ) {

    def result = 0

    chromosomeSizesFile.splitEachLine("\t") {fields ->

        def genomeId
        def path
	if ( fields[0].trim() == "${chromosome}" ) {
          //println "in if"
          result = fields[1].trim().toInteger()
          return 
        }

    }

    result
}


def getChromosomeChunkPairs ( chromosomeSize, chunkSize=5000000 ) {

  def result = []
  def numberOfChunks = chromosomeSize / chunkSize
  def remainder = chromosomeSize % chunkSize
 
  1.upto(numberOfChunks) {
    endPosition = it * chunkSize
    startPosition = (endPosition - chunkSize) + 1
    result = result + [[startPosition, endPosition ]]
  }

  if ( remainder > 0 ) {
    result = result + [[endPosition + 1 , endPosition + remainder ]]
  }
  
  result
}
