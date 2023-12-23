import os
from dotenv import load_dotenv
from pprint import pprint
import sys
load_dotenv(".env")

SECRET_KEY = 'ju=n4om3z00jd1+y2(ufn)g^@w-dj*&-45&4yd1_aiun50b6by'
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

sys.path.append(os.path.abspath(os.path.join(BASE_DIR,'ribctl')))       #! hack until ribctl is a separate pypi project
# sys.path.append(os.path.abspath(os.path.join(BASE_DIR,'ribctl','lib'))) 

        

# SECURITY WARNING: keep the secret key used in production secret!
# this should be either docker-mounted or populated through the utils module
RIBETL_DATA = os.environ["RIBETL_DATA"] if os.environ["RIBETL_DATA"] else os.path.join(BASE_DIR, "ribetldata")
# ※ Dont' forget to export pymol path ( we want to ship a built pymol, but python needs to be aware of it)
# export PYMOL_PATH=/home/rxz/dev/pymol3.11 && export PYTHONPATH="$PYTHONPATH:$PYMOL_PATH/modules/:"


vars          = ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD","NEO4J_CURRENTDB", "RIBETL_DATA"]
NEO4J_URI       :str= os.getenv("NEO4J_URI")
NEO4J_PASSWORD  :str= os.getenv("NEO4J_PASSWORD")
NEO4J_USER      :str= os.getenv("NEO4J_USER")
NEO4J_CURRENTDB :str= os.getenv("NEO4J_CURRENTDB")

for var in vars:
    print(var,":",os.getenv(var))
    if var not in os.environ:
        print("Environment variable {} not set".format(var))
        exit(1)


DEBUG = True
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
    'ninja',
    'mod_comp',
    'mod_db',
]
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


ROOT_URLCONF = 'rbxz_bend.urls'

TEMPLATES = [
    {
        'BACKEND': 'django.template.backends.django.DjangoTemplates',
        'DIRS': [ 'static/templates'],
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

WSGI_APPLICATION = 'rbxz_bend.wsgi.application'

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
TIME_ZONE = 'UTC'
USE_I18N = True
USE_L10N = True
USE_TZ = True



# def ready():
#     logs_dir = os.path.join(BASE_DIR, 'logs')
#     if not os.path.exists(logs_dir):
#         os.makedirs(logs_dir)

# ready()


NINJA_DOCS_VIEW  = 'swagger'
STATIC_URL       = '/static/'
STATIC_ROOT      = os.path.join(BASE_DIR, 'staticfiles')
STATICFILES_DIRS = [
    os.path.join(BASE_DIR, 'static'),
]
