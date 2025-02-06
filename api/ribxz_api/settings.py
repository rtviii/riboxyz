import os
from pathlib import Path
from dotenv import load_dotenv
import sys
load_dotenv(".env")

SECRET_KEY  = os.environ.get("SECRET_KEY")
DEBUG       = bool(os.environ.get("DEBUG", default=0))
BASE_DIR    = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
STATIC_URL = '/static/'
STATIC_ROOT = os.environ.get("DJANGO_STATIC_FILES_PATH", os.path.join(BASE_DIR, 'staticfiles'))
RIBCTL      = os.path.abspath(os.path.join(Path(BASE_DIR).parent.absolute()))

# RIBCTL     = os.path.abspath(os.path.join(Path(BASE_DIR).parent.absolute(),'ribctl'))
# print("Sourced RIBCTL: {}".format(RIBCTL))
sys.path.append(RIBCTL)       #! hack until ribctl is a separate pypi project

RIBETL_DATA = os.environ.get("RIBETL_DATA")
if RIBETL_DATA is None:
    raise Exception("$RIBETL_DATA is not set in the environment.")

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
CORS_ORIGIN_WHITELIST = ("https://ribosome.xyz", "https://*", "http://*")
SECURE_CROSS_ORIGIN_OPENER_POLICY = None




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


if DEBUG:
    STATICFILES_DIRS = [os.path.join(BASE_DIR, 'static')]
else:
    STATICFILES_DIRS = []