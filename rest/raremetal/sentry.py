import os
from raven.contrib.flask import Sentry
from raven.versioning import fetch_git_sha

def init_app(app):
  # Start logging errors
  try:
    sha = fetch_git_sha(os.path.join(app.root_path, "../../"))
  except:
    sha = "no-git"

  release = "raremetal-server@{}".format(sha)

  sentry_dsn = app.config.get("SENTRY_DSN")
  if sentry_dsn is not None:
    app.config["SENTRY_CONFIG"] = {
      "dsn": sentry_dsn,
      "release": release
    }

    if app.config["SENTRY_ENV"] is not None:
      app.config["SENTRY_CONFIG"]["environment"] = app.config["SENTRY_ENV"]

    # This attaches sentry to current_app.extensions['sentry']
    Sentry(app, register_signal=False, wrap_wsgi=False)
