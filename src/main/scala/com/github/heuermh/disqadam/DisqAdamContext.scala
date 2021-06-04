package com.github.heuermh.disqadam

import grizzled.slf4j.Logging

import htsjdk.samtools.ValidationStringency

import org.bdgenomics.adam.converters.{ AlignmentConverter, VariantContextConverter }

import org.bdgenomics.adam.ds.ADAMContext

import org.bdgenomics.adam.ds.read.AlignmentDataset

import org.bdgenomics.adam.ds.variant.VariantContextDataset

import org.disq_bio.disq.{
  FileCardinalityWriteOption,
  HtsjdkReadsRdd,
  HtsjdkReadsRddStorage,
  HtsjdkVariantsRdd,
  HtsjdkVariantsRddStorage,
  ReadsFormatWriteOption,
  VariantsFormatWriteOption
}

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
   * Load a path in BAM format as an AlignmentDataset with Disq.
   *
   * @param path path to load
   * @return the path in BAM format as an AlignmentDataset
   */
  def loadDisqBam(path: String): AlignmentDataset = {
    loadDisqBam(path, useNio = false, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in BAM format as an AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param useNio true to use nio, defaults to false
   * @return the path in BAM format as an AlignmentDataset
   */
  def loadDisqBam(path: String, useNio: Boolean): AlignmentDataset = {
    loadDisqBam(path, useNio, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in BAM format as an AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param useNio true to use nio, defaults to false
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   * @return the path in BAM format as an AlignmentDataset
   */
  def loadDisqBam(path: String, useNio: Boolean, stringency: ValidationStringency): AlignmentDataset = {
    logger.info(s"Loading $path in BAM format as AlignmentDataset with Disq...")
    logger.info(s"Using nio $useNio and stringency $stringency");

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
   * Save an AlignmentDataset to a path in BAM format with Disq.
   *
   * @param alignments AlignmentDataset to save
   * @param path path to save alignments to
   */
  def saveDisqBam(alignments: AlignmentDataset, path: String): Unit = {
    saveDisqBam(alignments, path, ValidationStringency.LENIENT)
  }

  /**
   * Save an AlignmentDataset to a path in BAM format with Disq.
   *
   * @param alignments AlignmentDataset to save
   * @param path path to save alignments to
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   */
  def saveDisqBam(alignments: AlignmentDataset, path: String, stringency: ValidationStringency) = {
    logger.info(s"Saving AlignmentDataset to $path in BAM format with Disq...")
    logger.info(s"Using stringency $stringency");

    val (header, records) = alignments.convertToSam()

    val htsjdkReadsRddStorage = HtsjdkReadsRddStorage
      .makeDefault(ac.sc)
      .validationStringency(stringency)

    val htsjdkReadsRdd = new HtsjdkReadsRdd(header, records.map(_.get()))

    htsjdkReadsRddStorage.write(htsjdkReadsRdd, path, ReadsFormatWriteOption.BAM, FileCardinalityWriteOption.SINGLE)
  }

  /**
   * Load a path in CRAM format as an AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param referencePath reference path
   * @return the path in BAM format as an AlignmentDataset
   */
  def loadDisqCram(path: String, referencePath: String): AlignmentDataset = {
    loadDisqCram(path, referencePath, useNio = false, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in CRAM format as an AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param referencePath reference path
   * @param useNio true to use nio, defaults to false
   * @return the path in BAM format as an AlignmentDataset
   */
  def loadDisqCram(path: String, referencePath: String, useNio: Boolean): AlignmentDataset = {
    loadDisqCram(path, referencePath, useNio, stringency = ValidationStringency.LENIENT)
  }

  /**
   * Load a path in CRAM format as an AlignmentDataset with Disq.
   *
   * @param path path to load
   * @param referencePath reference path
   * @param useNio true to use nio, defaults to false
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   * @return the path in CRAM format as an AlignmentDataset
   */
  def loadDisqCram(path: String, referencePath: String, useNio: Boolean, stringency: ValidationStringency): AlignmentDataset = {
    logger.info(s"Loading $path in CRAM format as AlignmentDataset with Disq...")
    logger.info(s"Using reference path $referencePath, nio $useNio, and stringency $stringency");

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
   * Save an AlignmentDataset to a path in CRAM format with Disq.
   *
   * @param alignments AlignmentDataset to save
   * @param path path to save alignments to
   * @param referencePath reference path
   */
  def saveDisqCram(alignments: AlignmentDataset, path: String, referencePath: String): Unit = {
    saveDisqCram(alignments, path, referencePath, ValidationStringency.LENIENT)
  }

  /**
   * Save an AlignmentDataset to a path in CRAM format with Disq.
   *
   * @param alignments AlignmentDataset to save
   * @param path path to save alignments to
   * @param referencePath reference path
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   */
  def saveDisqCram(alignments: AlignmentDataset, path: String, referencePath: String, stringency: ValidationStringency) = {
    logger.info(s"Saving AlignmentDataset to $path in CRAM format with Disq...")
    logger.info(s"Using reference path $referencePath and stringency $stringency");

    val (header, records) = alignments.convertToSam()

    val htsjdkReadsRddStorage = HtsjdkReadsRddStorage
      .makeDefault(ac.sc)
      .referenceSourcePath(referencePath)
      .validationStringency(stringency)

    val htsjdkReadsRdd = new HtsjdkReadsRdd(header, records.map(_.get()))

    htsjdkReadsRddStorage.write(htsjdkReadsRdd, path, ReadsFormatWriteOption.CRAM, FileCardinalityWriteOption.SINGLE)
  }

  /**
   * Load a path in VCF format as an VariantContextDataset with Disq.
   *
   * @param path path to load
   * @return the path in VCF format as an VariantContextDataset
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

  /**
   * Save a VariantContextDataset to a path in uncompressed VCF format with Disq.
   *
   * @param variantContexts Variant contexts to save
   * @param path path to save variants to
   */
  def saveDisqVcf(variantContexts: VariantContextDataset, path: String): Unit = {
    saveDisqVcf(variantContexts, path, ValidationStringency.LENIENT)
  }

  /**
   * Save a VariantContextDataset to a path in uncompressed VCF format with Disq.
   *
   * @param variantContexts Variant contexts to save
   * @param path path to save variants to
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   */
  def saveDisqVcf(variantContexts: VariantContextDataset, path: String, stringency: ValidationStringency) = {
    logger.info(s"Saving VariantContextDataset to $path in uncompressed VCF format with Disq...")
    logger.info(s"Using stringency $stringency");

    val (header, records) = variantContexts.convertToVcf(stringency)

    val htsjdkVariantsRddStorage = HtsjdkVariantsRddStorage.makeDefault(ac.sc)
    val htsjdkVariantsRdd = new HtsjdkVariantsRdd(header, records)

    htsjdkVariantsRddStorage.write(htsjdkVariantsRdd, path, VariantsFormatWriteOption.VCF, FileCardinalityWriteOption.SINGLE)
  }

  /**
   * Save a VariantContextDataset to a path in BGZF-compressed VCF format with Disq.
   *
   * @param variantContexts Variant contexts to save
   * @param path path to save variants to
   */
  def saveDisqVcfBgzf(variantContexts: VariantContextDataset, path: String): Unit = {
    saveDisqVcfBgzf(variantContexts, path, ValidationStringency.LENIENT)
  }

  /**
   * Save a VariantContextDataset to a path in BGZF-compressed VCF format with Disq.
   *
   * @param variantContexts Variant contexts to save
   * @param path path to save variants to
   * @param stringency validation stringency, defaults to ValidationStringency.LENIENT
   */
  def saveDisqVcfBgzf(variantContexts: VariantContextDataset, path: String, stringency: ValidationStringency) = {
    logger.info(s"Saving VariantContextDataset to $path in BGZF-compressed VCF format with Disq...")
    logger.info(s"Using stringency $stringency");

    val (header, records) = variantContexts.convertToVcf(stringency)

    val htsjdkVariantsRddStorage = HtsjdkVariantsRddStorage.makeDefault(ac.sc)
    val htsjdkVariantsRdd = new HtsjdkVariantsRdd(header, records)

    htsjdkVariantsRddStorage.write(htsjdkVariantsRdd, path, VariantsFormatWriteOption.VCF_BGZ, FileCardinalityWriteOption.SINGLE)
  }
}
