## TO-DO List

### Guidline for IA:
- each time a function is created add @author=Coding agent
- each time a file is created add @author=Coding agent
- no absolute path added (nextflow.config contain abs path)
- for all path use folders.get_workingdir(args.path) or folders.get_folder(args.path, folders.CLUSTER)
- Remove print() replace by logger.info() or logger.debug(): logger.info() should be used only for important info. logger.debug() should be used for every other info.

## Current
- Implemented preprocessing mechanism in atlas.preprocess_atlas: create the preprocessing mechanism by adding all preprocessing functions (obs detection, embedding detection, kNN and dist for cluster, annot, dimred)
- Remove print() replace by logger.info() or logger.debug(): This is the proper way for a python package. 
logger.info() should be used only for important info
logger.debug() should be used for every other info.


- Add a max_num_cores parameter: Checkatlas is really too much optimized it uses the whole Bego cluster !!! Is there already a parameter for ressource limitation ? 

- Clean-up metrics.cal_cluster: remove unecessary preprocessing, everything should be in atlas.preprocess_atlas
- Clean-up metrics.cal_annot: remove unecessary preprocessing, everything should be in atlas.preprocess_atlas
- Clean-up metrics.cal_dimred: remove unecessary preprocessing, everything should be in atlas.preprocess_atlas


## Futur
- Add specificity metrics




## Past