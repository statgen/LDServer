#!/usr/bin/env python3
import time
import requests
import psutil as ps
import os
from subprocess import Popen, check_call
from pathlib import Path

# Set number of workers
N_WORKERS = 1

# Set database location
# DB_NAME = "test_gunicorn.db"
# os.environ["SQLALCHEMY_DATABASE_URI"] = f"sqlite:///{DB_NAME}"

def has_n_children(pid, n):
  p = ps.Process(pid)
  kids = p.children()
  return len(kids) == n

# Need to call some flask command first due to an issue w/ sqlalchemy & create_all()
#check_call("flask show-references", shell=True)

# Start gunicorn
print("Starting gunicorn smoke test...")
proc = Popen(f"gunicorn --access-logfile - --error-logfile - -k gthread -w {N_WORKERS} --pythonpath rest 'ldserver:create_app()'", shell=True)

# Wait and check for workers to start
print("Waiting for workers...")
found_workers = False
guni_master_pid = ps.Process(proc.pid).children()[0].pid
wait = 1
for _ in range(5):
  if has_n_children(guni_master_pid, N_WORKERS):
    print(f"... found {N_WORKERS} workers")
    found_workers = True
    break
  else:
    time.sleep(wait)
    wait *= 2

if not found_workers:
  raise Exception("... failed to find workers")

# Test endpoint
print("Testing gunicorn endpoint...")
wait = 1
for _ in range(5):
  try:
    resp = requests.get("http://127.0.0.1:8000/correlations", timeout=(3, 3))
    if resp.ok:
      break
  except:
    time.sleep(wait)
    wait *= 2

if not resp.ok:
  raise Exception("During gunicorn smoke test: correlations endpoint did not return OK")

# Check for a simple key
if not 'data' in resp.json():
  raise Exception("During gunicorn smoke test: correlation endpoint did not return expected key 'data'")

# Terminate gunicorn
proc.kill()

# Remove test database
# Path("rest/ldserver/").joinpath(DB_NAME).unlink()

print("... endpoint tested successfully")
