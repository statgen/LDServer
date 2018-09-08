from flask import Flask

def create_app(test_config = None):
    app = Flask(__name__, instance_relative_config = True)

    if test_config is None:
        app.config.from_object('config.default')
        app.config.from_pyfile('config.py', silent = True)
        app.config.from_envvar('RESTLD_CONFIG_FILE', silent = True)
    else:
        app.config.from_mapping(test_config)

    from playground import web
    app.register_blueprint(web.bp)

    return app