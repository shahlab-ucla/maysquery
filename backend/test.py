import urllib.request
import json

data = json.dumps({"mz": 115.0037, "mode": "negative", "adducts": [], "tolerance_ppm": 5.0}).encode('utf-8')
req = urllib.request.Request('http://127.0.0.1:8002/api/run_pipeline', data=data, headers={'Content-Type': 'application/json'})
try:
    res = urllib.request.urlopen(req)
    print("SUCCESS:")
    print(res.read().decode())
except Exception as e:
    print("ERROR:")
    print(e)
