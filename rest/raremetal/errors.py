from __future__ import print_function
import traceback
import re
from flask import current_app, request, jsonify
from .sentry import sentry

class FlaskException(Exception):
  status_code = 400

  def __init__(self, message, status_code=None, payload=None):
    Exception.__init__(self)
    self.message = message
    if status_code is not None:
      self.status_code = status_code
    self.payload = payload

  def to_dict(self):
    rv = dict(self.payload or ())
    rv['message'] = self.message
    return rv

def handle_all(error):
  # Log all exceptions to Sentry, unless we're in a test environment like travis.
  if current_app.config["SENTRY_ENV"] != "travis":
    if sentry is not None:
      print("Attempting to log exception to Sentry...")
      sentry.captureException()
    else:
      print("Warning: Sentry not setup to log exception")

  # If we're in debug mode, re-raise the exception so we get the
  # browser debugger
  if current_app.debug:
    raise error

  # Also log the exception to the console.
  print("Exception thrown while handling request: " + request.url)
  traceback.print_exc()

  if isinstance(error, FlaskException):
    message = error.message
    code = error.status_code
  else:
    message = "An exception was thrown while handling the request. If you believe this request should have succeeded, please create an issue: https://github.com/statgen/LDServer/issues"
    code = 500

  # A little extra work to figure out the true request URL.
  # Requires the following set in apache:
  #   SetEnvIf Request_URI "^(.*)$" REQUEST_URI=$1
  #   RequestHeader set X-Request-Uri "%{REQUEST_URI}e"
  full_url = request.url
  real_uri = request.headers.get("X-Request-Uri")
  if real_uri is not None:
    match = re.search("\/(?P<api>\w+)\/(?P<version>v\d+)", real_uri)
    if match:
      api_name, api_version = match.groups()
      full_url = full_url.replace("/" + api_version,"/" + api_name + "/" + api_version)

  response = jsonify({
    "message": message,
    "request": full_url
  })
  response.status_code = code
  return response

def init_app(app):
  app.register_error_handler(Exception, handle_all)
