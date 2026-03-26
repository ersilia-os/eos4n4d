# Weighted Gram-Negative Accumulation Prediction

Prediction of small-molecule accumulation in Gram-negative bacteria using the weighted ensemble described by Köllen et al. The model combines an AutoML predictor trained on the Richter accumulation dataset with the TwinBooster zero-shot predictor, weighting both components by top-k Tanimoto similarity to the training set. It returns a single weighted accumulation score per compound for prioritization of molecules likely to accumulate in Gram-negative bacteria.

This model was incorporated on 2026-03-24.


## Information
### Identifiers
- **Ersilia Identifier:** `eos4n4d`
- **Slug:** `gram-negative-accumulation`

### Domain
- **Task:** `Annotation`
- **Subtask:** `Property calculation or prediction`
- **Biomedical Area:** `ADMET`
- **Target Organism:** `Escherichia coli`
- **Tags:** `E.coli`, `Fingerprint`, `Similarity`

### Input
- **Input:** `Compound`
- **Input Dimension:** `1`

### Output
- **Output Dimension:** `1`
- **Output Consistency:** `Fixed`
- **Interpretation:** Higher values indicate higher predicted Gram-negative accumulation.

Below are the **Output Columns** of the model:
| Name | Type | Direction | Description |
|------|------|-----------|-------------|
| weighted_accumulation | float | high | Weighted accumulation score |


### Source and Deployment
- **Source:** `Local`
- **Source Type:** `External`
- **S3 Storage**: [https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4n4d.zip](https://ersilia-models-zipped.s3.eu-central-1.amazonaws.com/eos4n4d.zip)

### Resource Consumption
- **Model Size (Mb):** `545`
- **Environment Size (Mb):** `8374`


### References
- **Source Code**: [https://github.com/sieber-lab/AIbiotics](https://github.com/sieber-lab/AIbiotics)
- **Publication**: [https://doi.org/10.1021/jacsau.5c00602](https://doi.org/10.1021/jacsau.5c00602)
- **Publication Type:** `Peer reviewed`
- **Publication Year:** `2025`
- **Ersilia Contributor:** [arnaucoma24](https://github.com/arnaucoma24)

### License
This package is licensed under a [GPL-3.0](https://github.com/ersilia-os/ersilia/blob/master/LICENSE) license. The model contained within this package is licensed under a [MIT](LICENSE) license.

**Notice**: Ersilia grants access to models _as is_, directly from the original authors, please refer to the original code repository and/or publication if you use the model in your research.


## Use
To use this model locally, you need to have the [Ersilia CLI](https://github.com/ersilia-os/ersilia) installed.
The model can be **fetched** using the following command:
```bash
# fetch model from the Ersilia Model Hub
ersilia fetch eos4n4d
```
Then, you can **serve**, **run** and **close** the model as follows:
```bash
# serve the model
ersilia serve eos4n4d
# generate an example file
ersilia example -n 3 -f my_input.csv
# run the model
ersilia run -i my_input.csv -o my_output.csv
# close the model
ersilia close
```

## About Ersilia
The [Ersilia Open Source Initiative](https://ersilia.io) is a tech non-profit organization fueling sustainable research in the Global South.
Please [cite](https://github.com/ersilia-os/ersilia/blob/master/CITATION.cff) the Ersilia Model Hub if you've found this model to be useful. Always [let us know](https://github.com/ersilia-os/ersilia/issues) if you experience any issues while trying to run it.
If you want to contribute to our mission, consider [donating](https://www.ersilia.io/donate) to Ersilia!
