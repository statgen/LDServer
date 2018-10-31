from flask import current_app, Blueprint, request, jsonify, make_response, abort, render_template

bp = Blueprint('api', __name__)

@bp.route('/')
def home():
    return render_template('home.html', data = { 'max_limit': current_app.config['API_MAX_PAGE_SIZE'], 'api_base_url': current_app.config['API_BASE_URL'] })