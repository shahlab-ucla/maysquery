import json
from main import run_single_pipeline
from models import MetaboliteInput
import asyncio

res = asyncio.run(run_single_pipeline(MetaboliteInput(query_type='ec', ec_number='4.2.1.2')))
print(json.dumps(res, default=str))
