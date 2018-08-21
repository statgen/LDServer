from flask import Flask

def create_app(test_config = None):
    app = Flask(__name__, instance_relative_config = True)

    if test_config is None:
        app.config.from_pyfile('config.py', silent = True)
        app.config.from_envvar('RESTLD_CONFIG_FILE', silent = True)
    else:
        app.config.from_mapping(test_config)

    from rest.model import db, load_references_command
    db.init_app(app)
    app.cli.add_command(load_references_command)

    from rest import api
    app.register_blueprint(api.bp)

    return app