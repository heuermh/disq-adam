package com.github.heuermh.disqadam

import org.bdgenomics.adam.ds.ADAMContext._

import org.bdgenomics.adam.util.ADAMFunSuite

class DisqAdamContextSuite extends ADAMFunSuite {

  override val appName: String = "disq-adam"
  override val properties: Map[String, String] = Map(
    "spark.serializer" -> "org.apache.spark.serializer.KryoSerializer",
    "spark.kryo.registrator" -> "org.bdgenomics.adam.serialization.ADAMKryoRegistrator,org.disq_bio.disq.serializer.DisqKryoRegistrator",
    "spark.kryo.referenceTracking" -> "true",
    "spark.kryo.registrationRequired" -> "true"
  )

  sparkTest("apply with ADAMContext") {
    DisqAdamContext(sc)
  }

  sparkTest("ctr with ADAMContext") {
    new DisqAdamContext(sc)
  }

  sparkTest("load BAM from Hadoop filesystem local path") {
    val dc = DisqAdamContext(sc)
    val path = testFile("NA12878.alignedHg38.duplicateMarked.baseRealigned.bam")

    val alignments = dc.loadDisqBam(path)
    assert(alignments.rdd.count() === 61614)
  }

  sparkTest("load CRAM from Hadoop filesystem local path") {
    val dc = DisqAdamContext(sc)
    val path = testFile("CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10m-10m100.cram")
    val referencePath = testFile("human_g1k_v37.20.21.fasta.gz")
    testFile("human_g1k_v37.20.21.fasta.gz.fai")
    testFile("human_g1k_v37.20.21.fasta.gz.gzi")

    val alignments = dc.loadDisqCram(path, referencePath)
    assert(alignments.rdd.count() === 654)
  }

  sparkTest("load VCF from Hadoop filesystem local path") {
    val dc = DisqAdamContext(sc)
    val path = testFile("CEUTrio.20.21.gatk3.4.g.vcf.bgz")
    testFile("CEUTrio.20.21.gatk3.4.g.vcf.idx")

    val variantContexts = dc.loadDisqVcf(path)
    assert(variantContexts.rdd.count() === 20037)
  }

  sparkTest("load BAM from nio local path") {
    val dc = DisqAdamContext(sc)
    val path = testFile("NA12878.alignedHg38.duplicateMarked.baseRealigned.bam")

    val alignments = dc.loadDisqBam(path, useNio = true)
    assert(alignments.rdd.count() === 61614)
  }

  sparkTest("load CRAM from nio local path") {
    val dc = DisqAdamContext(sc)
    val path = testFile("CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10m-10m100.cram")
    val referencePath = testFile("human_g1k_v37.20.21.fasta.gz")
    testFile("human_g1k_v37.20.21.fasta.gz.fai")
    testFile("human_g1k_v37.20.21.fasta.gz.gzi")

    val alignments = dc.loadDisqCram(path, referencePath, useNio = true)
    assert(alignments.rdd.count() === 654)
  }
}
