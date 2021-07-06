from flask import Flask

def create_app(test_config = None):
    app = Flask(__name__, instance_relative_config = True)

    if test_config is None:
        app.config.from_object('config.default')
        app.config.from_pyfile('config.py', silent = True)
        app.config.from_envvar('RESTLD_CONFIG_FILE', silent = True)
    else:
        app.config.from_mapping(test_config)

    from ldserver.model import db, load_correlations, load_references_command, show_references_command, \
        add_reference_command, create_subset_command, show_genotypes_command, show_samples_command
    db.init_app(app)
    app.cli.add_command(load_references_command)
    app.cli.add_command(show_references_command)
    app.cli.add_command(add_reference_command)
    app.cli.add_command(create_subset_command)
    app.cli.add_command(show_genotypes_command)
    app.cli.add_command(show_samples_command)

    from ldserver import api
    app.register_blueprint(api.bp)

    if app.config['GZIP_COMPRESSION']:
        app.config['COMPRESS_MIMETYPES'] = ['application/json', 'application/msgpack']
        app.config['COMPRESS_LEVEL'] = 3
        app.config['COMPRESS_MIN_SIZE'] = 500
        api.compress.init_app(app)

    with app.app_context():
        app.config['CORRELATIONS'] = load_correlations()

    return app
