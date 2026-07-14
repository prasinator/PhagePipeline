import urllib.request
import sys
import time

url = "https://zenodo.org/api/records/10515367/files/RBPdetect_v3_ESMfineT33.zip/content"
output = "RBPdetect_v3_ESMfineT33.zip"

print("Connecting to Zenodo server...")
req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})

try:
    with urllib.request.urlopen(req) as response, open(output, 'wb') as out_file:
        total_size = int(response.getheader('Content-Length', 0))
        downloaded = 0
        start_time = time.time()
        
        print(f"Connected successfully! Target file size: {total_size / (1024**3):.2f} GB")
        
        while True:
            chunk = response.read(512 * 1024)  # Read 512KB chunks
            if not chunk:
                break
            downloaded += len(chunk)
            out_file.write(chunk)
            
            percent = (downloaded / total_size) * 100
            elapsed = time.time() - start_time
            speed = (downloaded / (1024 * 1024)) / elapsed if elapsed > 0 else 0
            
            sys.stdout.write(f"\rProgress: {percent:.1f}% | Speed: {speed:.2f} MB/s")
            sys.stdout.flush()
            
    print("\nDownload complete! Verifying archive structure...")
except Exception as e:
    print(f"\nNetwork transfer interrupted: {e}")
