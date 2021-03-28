package com.github.heuermh.disqadam

import grizzled.slf4j.Logging

import htsjdk.samtools.ValidationStringency

import org.bdgenomics.adam.converters.{ AlignmentConverter, VariantContextConverter }

import org.bdgenomics.adam.ds.ADAMContext

import org.bdgenomics.adam.ds.read.AlignmentDataset

import org.bdgenomics.adam.ds.variant.VariantContextDataset

import org.disq_bio.disq.{ HtsjdkReadsRddStorage, HtsjdkVariantsRddStorage }

object DisqAdamContext {

  /**
   * Create a new DisqAdamContext extending the specified ADAMContext.
   *
   * @param ac ADAMContext to extend.
   */
  def apply(ac: ADAMContext): DisqAdamContext = {
    new DisqAdamContext(ac)
  }
}

class DisqAdamContext(@transient val ac: ADAMContext) extends Serializable with Logging {

  /**
   * Load a path in BAM format as a AlignmentDataset with Disq.
   *
   * @param path path to load
   * @return the path in BAM format as a AlignmentDataset
   */
  def loadDisqBam(path: String): AlignmentDataset = {
    loadDisqBam(path, useNio = false, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in BAM format as a AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param useNio true to use nio, defaults to false
   * @return the path in BAM format as a AlignmentDataset
   */
  def loadDisqBam(path: String, useNio: Boolean): AlignmentDataset = {
    loadDisqBam(path, useNio, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in BAM format as a AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param useNio true to use nio, defaults to false
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   * @return the path in BAM format as a AlignmentDataset
   */
  def loadDisqBam(path: String, useNio: Boolean, stringency: ValidationStringency): AlignmentDataset = {
    logger.info(s"Loading $path in BAM format as AlignmentDataset with Disq...")
    val htsjdkReadsRddStorage = HtsjdkReadsRddStorage
      .makeDefault(ac.sc)
      .useNio(useNio)
      .validationStringency(stringency)
    val htsjdkReadsRdd = htsjdkReadsRddStorage.read(path)

    val header = htsjdkReadsRdd.getHeader
    val references = AlignmentConverter.references(header)
    val readGroups = AlignmentConverter.readGroups(header)
    val processingSteps = AlignmentConverter.processingSteps(header)
    val converter = new AlignmentConverter()

    val reads = htsjdkReadsRdd.getReads
    val alignmentRdd = reads.rdd.map(converter.convert)

    AlignmentDataset(alignmentRdd, references, readGroups, processingSteps)
  }

  /**
   * Load a path in BAM format as a AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param referencePath reference path
   * @return the path in BAM format as a AlignmentDataset
   */
  def loadDisqCram(path: String, referencePath: String): AlignmentDataset = {
    loadDisqCram(path, referencePath, useNio = false, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in BAM format as a AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param referencePath reference path
   * @param useNio true to use nio, defaults to false
   * @return the path in BAM format as a AlignmentDataset
   */
  def loadDisqCram(path: String, referencePath: String, useNio: Boolean): AlignmentDataset = {
    loadDisqCram(path, referencePath, useNio, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in CRAM format as a AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param referencePath reference path
   * @param useNio true to use nio, defaults to false
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   * @return the path in CRAM format as a AlignmentDataset
   */
  def loadDisqCram(path: String, referencePath: String, useNio: Boolean, stringency: ValidationStringency): AlignmentDataset = {
    logger.info(s"Loading $path in CRAM format as AlignmentDataset with Disq...")
    val htsjdkReadsRddStorage = HtsjdkReadsRddStorage
      .makeDefault(ac.sc)
      .referenceSourcePath(referencePath)
      .useNio(useNio)
      .validationStringency(stringency)
    val htsjdkReadsRdd = htsjdkReadsRddStorage.read(path)

    val header = htsjdkReadsRdd.getHeader
    val references = AlignmentConverter.references(header)
    val readGroups = AlignmentConverter.readGroups(header)
    val processingSteps = AlignmentConverter.processingSteps(header)
    val converter = new AlignmentConverter()

    val reads = htsjdkReadsRdd.getReads
    val alignmentRdd = reads.rdd.map(converter.convert)

    AlignmentDataset(alignmentRdd, references, readGroups, processingSteps)
  }

  /**
   * Load a path in VCF format as a VariantContextDataset with Disq.
   *
   * @param path path to load
   * @return the path in VCF format as a VariantContextDataset
   */
  def loadDisqVcf(path: String): VariantContextDataset = {
    loadDisqVcf(path, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in VCF format as a VariantContextDataset with Disq.
   *
   * @param path path to load
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   * @return the path in VCF format as a VariantContextDataset
   */
  def loadDisqVcf(path: String, stringency: ValidationStringency): VariantContextDataset = {
    logger.info(s"Loading $path in VCF format as VariantContextDataset with Disq...")
    val htsjdkVariantsRddStorage = HtsjdkVariantsRddStorage.makeDefault(ac.sc)
    val htsjdkVariantsRdd = htsjdkVariantsRddStorage.read(path)

    val header = htsjdkVariantsRdd.getHeader
    val headerLines = VariantContextConverter.headerLines(header)
    val samples = VariantContextConverter.samples(header)
    val references = VariantContextConverter.references(header)
    val converter = VariantContextConverter(headerLines, stringency, ac.sc.hadoopConfiguration)

    val variants = htsjdkVariantsRdd.getVariants
    val variantContextRdd = variants.rdd.flatMap(converter.convert)

    VariantContextDataset(
      variantContextRdd,
      references,
      samples,
      VariantContextConverter.cleanAndMixInSupportedLines(headerLines, stringency, logger.logger)
    )
  }
}
