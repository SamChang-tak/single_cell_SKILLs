# CellRank analysis report

## Objective

{{BIOLOGICAL_QUESTION}}

## Input data

- Dataset: {{DATASET}}
- Cells × genes: {{DIMENSIONS}}
- Cluster key: {{CLUSTER_KEY}}
- Directional signal and key: {{DIRECTION_SIGNAL}}
- Embedding: {{EMBEDDING}}

## Methods

- Preprocessing: {{PREPROCESSING}}
- Kernel(s): {{KERNELS}}
- Kernel weights and rationale: {{KERNEL_WEIGHTS}}
- Estimator and parameters: {{ESTIMATOR}}
- Initial/terminal-state definitions: {{STATE_DEFINITIONS}}
- Software versions: {{VERSIONS}}

## Quality checks

- Neighbor graph: {{NEIGHBOR_GRAPH_QC}}
- Transition matrix: {{TRANSITION_MATRIX_QC}}
- State stability/sensitivity: {{STATE_SENSITIVITY}}
- Fate-probability validation: {{FATE_QC}}

## Results

### Macrostates and terminal states

{{MACROSTATE_RESULTS}}

### Fate probabilities

{{FATE_RESULTS}}

### Lineage drivers and expression trends

{{DRIVER_AND_TREND_RESULTS}}

## Limitations

{{LIMITATIONS}}

## Output files

{{OUTPUT_MANIFEST}}
