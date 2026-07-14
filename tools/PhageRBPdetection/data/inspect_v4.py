import urllib.request
import json
import sys

url = "https://zenodo.org/api/records/14810759"
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64)'})

print("Fetching metadata from Zenodo record 14810759...")
try:
    with urllib.request.urlopen(req) as response:
        data = json.loads(response.read().decode())
    
    print("Record files:")
    for f in data.get('files', []):
        print(json.dumps(f, indent=2))

except Exception as e:
    print(f"An error occurred: {e}")
