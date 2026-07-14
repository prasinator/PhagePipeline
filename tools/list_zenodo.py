import urllib.request
import json

try:
    url = "https://zenodo.org/api/records/10515367"
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    response = urllib.request.urlopen(req)
    data = json.loads(response.read().decode())
    for f in data.get('files', []):
        print(f"File: {f.get('key')}")
        print(f"Link: {f.get('links', {}).get('self')}")
        print(f"Size: {f.get('size')} bytes")
        print("---")
except Exception as e:
    print(f"Error: {e}")
