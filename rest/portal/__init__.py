from flask import Flask

def create_app(test_config = None):
    app = Flask(__name__, instance_relative_config = True)

    if test_config is None:
        app.config.from_object('config.default')
        app.config.from_pyfile('config.py', silent = True)
        app.config.from_envvar('RESTLD_CONFIG_FILE', silent = True)
    else:
        app.config.from_mapping(test_config)

    from .model import db, load_correlations, load_genotype_datasets_command, show_genotype_datasets_command, \
        add_genotype_dataset_command, create_subset_command, show_genotypes_command, show_samples_command, \
        add_phenotypes_command
    
    db.init_app(app)
    app.cli.add_command(load_genotype_datasets_command)
    app.cli.add_command(show_genotype_datasets_command)
    app.cli.add_command(add_genotype_dataset_command)
    app.cli.add_command(add_phenotypes_command)
    app.cli.add_command(create_subset_command)
    app.cli.add_command(show_genotypes_command)
    app.cli.add_command(show_samples_command)

    from rest import api
    app.register_blueprint(api.bp)

    if app.config['GZIP_COMPRESSION']:
        app.config['COMPRESS_MIMETYPES'] = ['application/json']
        app.config['COMPRESS_LEVEL'] = 3
        app.config['COMPRESS_MIN_SIZE'] = 500
        api.compress.init_app(app)

    with app.app_context():
        app.config['CORRELATIONS'] = load_correlations()

    return app