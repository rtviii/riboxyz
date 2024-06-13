import os
from pathlib import Path
from dotenv import load_dotenv
import sys
load_dotenv(".env")

SECRET_KEY = os.environ.get("SECRET_KEY")
DEBUG      = bool(os.environ.get("DEBUG", default=0))
BASE_DIR   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RIBCTL     = os.path.abspath(os.path.join(Path(BASE_DIR).parent.absolute()))
# RIBCTL     = os.path.abspath(os.path.join(Path(BASE_DIR).parent.absolute(),'ribctl'))
# print("Sourced RIBCTL: {}".format(RIBCTL))
sys.path.append(RIBCTL)       #! hack until ribctl is a separate pypi project

RIBETL_DATA = os.environ["RIBETL_DATA"] if os.environ["RIBETL_DATA"] else os.path.join(BASE_DIR, "ribetldata")

# ※ Dont' forget to export pymol path ( we want to ship a built pymol, but python needs to be aware of it)
# export PYMOL_PATH=/home/rxz/dev/pymol3.11 && export PYTHONPATH="$PYTHONPATH:$PYMOL_PATH/modules/:"
ALLOWED_HOSTS = ["*"]
INSTALLED_APPS = [
    'django.contrib.admin',
    'django.contrib.auth',
    'django.contrib.contenttypes',
    'django.contrib.sessions',
    'django.contrib.messages',
    'django.contrib.staticfiles',
    'rest_framework',
    'corsheaders',
    'ninja']

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
    'corsheaders.middleware.CorsMiddleware',
    'django.middleware.common.CommonMiddleware'
]
CORS_ALLOW_ALL_ORIGINS = True # If this is used then `CORS_ALLOWED_ORIGINS` will not have any effect
CORS_ALLOW_CREDENTIALS = True
# CORS_ALLOWED_ORIGINS = [
#     'https://*',
#     'http://*',
#     'localhost',
# ]

ROOT_URLCONF = 'ribxz_api.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [ os.path.join(BASE_DIR,'templates')],
        'APP_DIRS': True,
        'OPTIONS': {
            'context_processors': [
                'django.template.context_processors.debug',
                'django.template.context_processors.request',
                'django.contrib.auth.context_processors.auth',
                'django.contrib.messages.context_processors.messages',
            ],
        },
    },
]

WSGI_APPLICATION = 'ribxz_api.wsgi.application'

DATABASES = {}

# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator', },
    {'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator', },
    {'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator', },
    {'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator', },
]

LANGUAGE_CODE = 'en-us'
TIME_ZONE     = 'UTC'
USE_I18N      = True
USE_L10N      = True
USE_TZ        = True



# NINJA_DOCS_VIEW  = 'swagger'

# STATIC_URL       = '/static/'
# STATIC_ROOT      = os.path.join(BASE_DIR, 'staticfiles')
# STATICFILES_DIRS = [ os.path.join(BASE_DIR, 'static'), ]

# STATIC_URL lets you namespace your static files to avoid url conflicts and make them inaccessible from the browser
STATIC_URL       = '/static/'
# STATIC_ROOT is where all the static files are collected by manage.py collectstatic
STATIC_ROOT      = '/srv/www/static'
STATICFILES_DIRS = [ os.path.join(BASE_DIR, 'static'), ]
