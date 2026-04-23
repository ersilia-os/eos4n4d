---
language:
- en
metrics:
- perplexity
pipeline_tag: sentence-similarity
author: Maximilian G. Schuh
tags:
- PubChem
- chemistry
- biology
- deberta-v3
---
# TwinBooster

## PubChemDeBERTa-augmented: Fine-tuned DeBERTa V3 base on PubChem bioassay corpora

[![arXiv](https://img.shields.io/badge/arXiv-2401.04478-b31b1b.svg)](https://arxiv.org/abs/2401.04478)

### Synergising Large Language Models with Barlow Twins and Gradient Boosting for Enhanced Molecular Property Prediction

Maximilian G. Schuh, Davide Boldini, Stephan A. Sieber

@ Chair of Organic Chemistry II,
TUM School of Natural Sciences,
Technical University of Munich

**Abstract**

The success of drug discovery and development relies on the precise prediction of molecular activities and properties. While in silico molecular property prediction has shown remarkable potential, its use has been limited so far to assays for which large amounts of data are available. In this study, we use a fine-tuned large language model to integrate biological assays based on their textual information, coupled with Barlow Twins, a Siamese neural network using a novel self-supervised learning approach. This architecture uses both assay information and molecular fingerprints to extract the true molecular information. TwinBooster enables the prediction of properties of unseen bioassays and molecules by providing state-of-the-art zero-shot learning tasks. Remarkably, our artificial intelligence pipeline shows excellent performance on the FS-Mol benchmark. This breakthrough demonstrates the application of deep learning to critical property prediction tasks where data is typically scarce. By accelerating the early identification of active molecules in drug discovery and development, this method has the potential to help streamline the identification of novel therapeutics. 
