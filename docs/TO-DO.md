## TO-DO List

### Guidline for IA:
- each time a function is created add @author=Coding agent
- each time a file is created add @author=Coding agent
- no absolute path added (nextflow.config contain abs path)
- for all path use folders.get_workingdir(args.path) or folders.get_folder(args.path, folders.CLUSTER)
- Remove print() replace by logger.info() or logger.debug(): logger.info() should be used only for important info. logger.debug() should be used for every other info.

## Current

### Checkatlas:
- Implemented preprocessing mechanism in atlas.preprocess_atlas: create the preprocessing mechanism by adding all preprocessing functions (obs detection, embedding detection, kNN and dist for cluster, annot, dimred) [DONE]
==> Implemented the preprocess_atlas function as @atlas.py

- Remove print() replace by logger.info() or logger.debug(): This is the proper way for a python package. [DONE]
logger.info() should be used only for important info
logger.debug() should be used for every other info.
==> Implmented this 

- Add a max_num_cores parameter: Checkatlas is really too much optimized it uses the whole Bego cluster !!! Is there already a parameter for ressource limitation ? [DONE]
==> Yeas, we are taking maximum 48 cores and also 80 percent of the idle GPU.

- Clean-up metrics.cal_cluster: remove unecessary preprocessing, everything should be in atlas.preprocess_atlas [DONE]
- Clean-up metrics.cal_annot: remove unecessary preprocessing, everything should be in atlas.preprocess_atlas
- Clean-up metrics.cal_dimred: remove unecessary preprocessing, everything should be in atlas.preprocess_atlas
==> Now these are cleaned 

- Clean code, format code, lint code : make fmt and make lint [NOT_DONE]
poetry install --with dev [HAVE TO DO ...]



### Nextflow:


## Futur
- Add specificity metrics
- Create test functions : make test and tests/ folder



## Past