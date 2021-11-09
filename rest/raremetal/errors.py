import traceback
import re
import werkzeug
import sys
from flask import current_app, request, jsonify
from core.pywrapper import LDServerGenericException

class FlaskException(Exception):
  status_code = 400

  def __init__(self, message, status_code=None, secret=None):
    """
    Construct an exception representing an error that occurred when running flask or a flask subcommand.
    This exception can include "secret" data that will not be included in a HTTP response message, but will be
    delivered to the server log.

    :param message: Message describing the exception. This will logged in the console, in Sentry, and in the HTTP response (if applicable.)
    :param status_code: HTTP status code.
    :param secret: Message with additional information, such as filepaths, that will only appear in the server log.
    """

    super().__init__(message)
    self.message = message
    self._secret = secret

    if status_code is not None:
      self.status_code = status_code

  @property
  def secret(self):
    return self._secret

  @secret.setter
  def secret(self, value):
    self._secret = value

def handle_all(error):
  sentry = current_app.extensions.get("sentry")

  # Try to log exception to Sentry if it is configured.
  if sentry is not None:
    print("Attempting to log exception to Sentry...")
    sentry.captureException()
  else:
    print("Sentry not setup to log exceptions")

  # If we're in debug mode, re-raise the exception so we get the
  # browser debugger
  if current_app.debug:
    raise error

  # Also log the exception to the console.
  print("Exception thrown while handling request: " + request.url, file=sys.stderr)
  traceback.print_exc() # defaults to stderr
  print(str(error), file=sys.stderr)
  if isinstance(error, FlaskException) and error.secret is not None:
    print(error.secret, file=sys.stderr)

  if isinstance(error, FlaskException):
    message = error.message
    code = error.status_code
  elif isinstance(error, werkzeug.exceptions.NotFound):
    message = error.description
    code = error.code
  elif hasattr(error, 'args') and (error.args is not None) and (len(error.args) > 0) and isinstance(error.args[0], LDServerGenericException):
    # This is a safe exception type with sensitive information kept in a separate string object.

    # Print extra C++ exception information to stderr
    print("C++ exception thrown:\n∟ public msg:  {}\n∟ private msg: {}".format(str(error), error.args[0].get_secret()), file=sys.stderr)

    # Raise a flask exception to tell the developer what went wrong
    message = str(error)
    code = 400 # Any C++ exception should be a HTTP 400
  else:
    message = "An unexpected error occurred on the server while processing the request. Please ask the admin to check the server logs for more information."
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
    "data": {},
    "error": message,
    "request": full_url
  })
  response.status_code = code
  return response

def init_app(app):
  app.register_error_handler(Exception, handle_all)
