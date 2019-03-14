from invoke import task

"""
LDServer tasks

To run:

  * invoke build
  * invoke test
  * invoke version --part <major|minor|part|beta|dev>
"""

VALID_PARTS = "major minor patch beta dev".split()

@task
def build(ctx):
  """
  Recompile the C++ component (pywrapper)
  """

  ctx.run("cget install --update core")

@task(build)
def test(ctx):
  """
  Run all test cases using tox. This task will automatically execute build.
  """

  ctx.run("tox")

@task(test)
def version(ctx, part):
  """
  Bump the version number using semantic versioning scheme. The configuration is stored in setup.cfg.

  The version will only be bumped if the build and test tasks run successfully.
  """

  if part not in VALID_PARTS:
    raise ValueError("Invalid part specifier: {}, must be one of {}".format(part, ", ".join(VALID_PARTS)))

  ctx.run("bumpversion --verbose {}".format(part))
