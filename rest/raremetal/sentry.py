import os
from raven.contrib.flask import Sentry
from raven.versioning import fetch_git_sha

sentry = None

def init_app(app):
  global sentry

  # Start logging errors
  sha = fetch_git_sha(os.path.join(app.root_path, "../../"))
  release = "raremetal-server@{}".format(sha)

  if sentry is None:
    if "SENTRY_DSN" in app.config:
      app.config["SENTRY_CONFIG"] = {
        "dsn": app.config["SENTRY_DSN"],
        "release": release
      }

      if app.config["SENTRY_ENV"] is not None:
        app.config["SENTRY_CONFIG"]["environment"] = app.config["SENTRY_ENV"]

      sentry = Sentry(app, register_signal=False, wrap_wsgi=False)