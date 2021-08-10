#!/usr/bin/env python3
from sarge import run, Capture
import time
import requests

# Start gunicorn
with Capture(buffer_size=-1) as c:
  with run("gunicorn --access-logfile - --error-logfile - -k gthread --pythonpath rest 'ldserver:create_app()'", async_=True, stdout=c, stderr=c) as p:
    # Check to make sure a worker was booted
    print("Waiting for worker to boot...")
    m = c.expect("Booting worker with pid", timeout=3)
    if m is None:
      raise Exception("During gunicorn smoke test: never saw expected gunicorn worker booting")
    else:
      print("... found worker")

    # Test endpoint
    print("Testing gunicorn endpoint...")
    resp = requests.get("http://127.0.0.1:8000/correlations", timeout=(3, 3))
    if not resp.ok:
      raise Exception("During gunicorn smoke test: correlations endpoint did not return OK")

    # Check for a simple key
    if not 'data' in resp.json():
      raise Exception("During gunicorn smoke test: correlation endpoint did not return expected key 'data'")

    print("... endpoint tested successfully")

    # We're all set, can shut down gunicorn
    p.commands[0].terminate()
