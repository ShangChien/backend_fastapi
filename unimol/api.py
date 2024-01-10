from numpy import ndarray
from .utils import gen3D,dataUnimol,RES
from typing import Any
from unimol_tools import MolPredict, UniMolRepr
from pathlib import Path as P
import logging
import copy

logger = logging.getLogger(__name__)

models={
    'homo': MolPredict(load_model=P(P(__file__).parent, './jobs/homo/exp')),
    'lumo': MolPredict(load_model=P(P(__file__).parent, './jobs/lumo/exp')),
    'eg'  : MolPredict(load_model=P(P(__file__).parent, './jobs/eg/exp')),
    't1'  : MolPredict(load_model=P(P(__file__).parent, './jobs/t1/exp')),
    's1'  : MolPredict(load_model=P(P(__file__).parent, './jobs/s1/exp')),
    'repr': UniMolRepr(data_type='molecule')
}
    
def inference(data: dataUnimol) -> RES[dataUnimol]:
    global models
    logger.info('input data:', data.model_dump_json())
    inputData = {}
    RES_result = RES[dataUnimol]()
    RES_result.data = copy.deepcopy(data)

    # 前处理获取3d数据
    if data.atoms and data.coordinates:
        logger.info('process atoms and coord')
        inputData['atoms'] = data.atoms
        inputData['coordinates'] = data.coordinates
    elif data.molBlocks or data.smiles:
        logger.info('before generate 3D by RDKit')
        res_3D:RES[dataUnimol] = gen3D(molBlocks=data.molBlocks, smiles=data.smiles)
        if res_3D.success and res_3D.data is not None:
            inputData['atoms'] = res_3D.data.atoms
            inputData['coordinates'] = res_3D.data.coordinates
            logger.info('success to get 3D info')
        else:
            error_str:str = f'fail to get 3D info before unimol: {res_3D.error}'
            logger.error(error_str)
            RES_result.error = RES_result.error + error_str
            RES_result.success = False
            return RES_result
    else:
        error_str:str = 'input data wrong: data is invalid'
        logger.error(error_str)
        RES_result.error = RES_result.error + error_str
        RES_result.success = False
        return RES_result
    logger.info("get 3D data", data.model_dump_json())

    # unimol推理
    if data.models:
        for model in data.models:
            model:str = model.lower()
            results:list[Any] = []
            if model == 'repr':
                reprs:dict = models['repr'].get_repr(inputData)
                results = reprs["cls_repr"]
            elif model in models.keys():
                _results:ndarray = models[model].predict(data=inputData)
                results = _results.flatten().tolist()
            else:
                logger.error(f'unimol model type error: {model}')
            RES_result.data.results[model] = results
            logger.info(f'success run unimol inference: {model}; total mol-num: {len(results)}')
    else:
        error_str:str = 'input data wrong: models is invalid'
        logger.error(error_str)
        RES_result.error = RES_result.error + error_str
        RES_result.success = False
    return RES_result