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

def value_ok(v):
  if hasattr(v, "__len__"):
    if len(v) > 0:
      return v
  elif v is not None:
    return v

def retry(func, *args, **kwargs):
  init_wait = 1
  multiply = 2
  tries = 5

  wait = init_wait
  for _ in range(tries):
    try:
      v = func(*args, **kwargs)
      if value_ok(v):
        return v
      else:
        time.sleep(wait)
        wait *= multiply
    except:
      time.sleep(wait)
      wait *= multiply

  raise Exception(f"func {str(func)} failed after {tries} attempts")

def has_n_children(pid, n):
  print(f"Expecting {n} kids for {pid}")
  p = ps.Process(pid)
  kids = p.children()
  return len(kids) == n

def find_guni_pid(pid):
  return ps.Process(pid).children()[0].pid

# Need to call some flask command first due to an issue w/ sqlalchemy & create_all()
#check_call("flask show-references", shell=True)

# Start gunicorn
print("Starting gunicorn smoke test...")
proc = Popen(f"gunicorn --access-logfile - --error-logfile - -k gthread -w {N_WORKERS} --pythonpath rest 'ldserver:create_app()'", shell=True)

# Find gunicorn master PID
print(f"Main PID: {proc.pid}")
guni_master_pid = retry(find_guni_pid, proc.pid)
print(f"Master gunicorn PID: {guni_master_pid}")

# Wait and check for workers to start
print("Waiting for workers...")
found_workers = retry(has_n_children, guni_master_pid, N_WORKERS)

if not found_workers:
  raise Exception("... failed to find workers")

# Test endpoint
print("Testing gunicorn endpoint...")
resp = retry(requests.get, "http://127.0.0.1:8000/correlations", timeout=(3, 3))

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
