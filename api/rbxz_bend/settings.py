import datetime
import os
import logging.handlers
import logging
import typing

class DatedRotatingFileHandler(logging.handlers.RotatingFileHandler):
    """
    Handler for rotating log files with a date-based filename.
    """
    def __init__(self, filename, mode='a', maxBytes=0, backupCount=0, encoding=None, delay=False):
        # Append the current datetime to the log filename
        now = datetime.datetime.now().strftime('%Y-%m-%d_%H-%M-%S')
        filename = os.path.splitext(filename)[0] + '_' + now + os.path.splitext(filename)[1] + '.log'
        super().__init__(filename, mode, maxBytes, backupCount, encoding, delay)

LOGGERS = typing.Literal['general', 'updates', 'accesses','computations']
def get_logger(loggername: LOGGERS )-> logging.Logger:
    """if not one of the default loggers is specified, then a rotating log with a new date is created"""
    if loggername not in typing.get_args(LOGGERS):
        # Create a custom logger and handler instance
        logger = logging.getLogger(loggername)
        handler = DatedRotatingFileHandler(os.path.join(BASE_DIR,'logs',loggername), maxBytes=1024*1024*5, backupCount=5)
        handler.setLevel(logging.DEBUG)
        formatter = logging.Formatter('%(asctime)s - %(filename)s - %(levelname)s - %(message)s')
        handler.setFormatter(formatter)
        logger.addHandler(handler)
        logger.setLevel(logging.DEBUG)
        return logger
    else:
        return logging.getLogger(loggername)
        



# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'ju=n4om3z00jd1+y2(ufn)g^@w-dj*&-45&4yd1_aiun50b6by'
BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# this should be either docker-mounted or populated through the utils module
RIBETL_DATA = os.environ["RIBETL_DATA"] if os.environ["RIBETL_DATA"] else os.path.join(
    BASE_DIR, "ribetldata")

# ※ Dont' forget to export pymol path ( we want to ship a built pymol, but python needs to be aware of it)
# export PYMOL_PATH=/home/rxz/dev/pymol3.11 && export PYTHONPATH="$PYTHONPATH:$PYMOL_PATH/modules/:"

vars = ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD",
        "NEO4J_CURRENTDB", "RIBETL_DATA"]
print(" ---------------------------------------------- App has been reset. -----------------------------------------------")
for var in vars:
    # print("Environment variable {}:\t{}".format( var, os.getenv(var) ))
    if var not in os.environ:
        print("Environment variable {} not set".format(var))
        exit(1)

NEO4J_URI = os.getenv("NEO4J_URI")
NEO4J_PASSWORD = os.getenv("NEO4J_PASSWORD")
NEO4J_USER = os.getenv("NEO4J_USER")
NEO4J_CURRENTDB = os.getenv("NEO4J_CURRENTDB")

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
    'ninja',
    'mod_comp',
    'mod_db',
    'mod_utils',

]
NINJA_DOCS_VIEW='swagger'

MIDDLEWARE = [
    'django.middleware.security.SecurityMiddleware',
    'django.contrib.sessions.middleware.SessionMiddleware',
    'django.middleware.common.CommonMiddleware',
    'django.middleware.csrf.CsrfViewMiddleware',
    'django.contrib.auth.middleware.AuthenticationMiddleware',
    'django.contrib.messages.middleware.MessageMiddleware',
    'django.middleware.clickjacking.XFrameOptionsMiddleware',
]

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



def ready():
    logs_dir = os.path.join(BASE_DIR, 'logs')
    if not os.path.exists(logs_dir):
        os.makedirs(logs_dir)

ready()

LOGGING = {
    'version': 1,
    'disable_existing_loggers': False,
    'formatters': {
        'standard': {
            'format': '%(asctime)s - %(filename)s - %(levelname)s - %(message)s',
        },
    },
    'handlers': {
        'general': {
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': os.path.join(BASE_DIR,'logs','general.log'),
            'maxBytes': 1024 * 1024 * 5, # 5MB
            'backupCount': 5,
            'formatter': 'standard',
        },
        'updates': {
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': os.path.join(BASE_DIR,'logs','updates.log'),
            'maxBytes': 1024 * 1024 * 5, # 5MB
            'backupCount': 5,
            'formatter': 'standard',
        },
        'accesses': {
            'class'      : 'logging.handlers.RotatingFileHandler',
            'filename'   : os.path.join(BASE_DIR,'logs','accesses.log'),
            'maxBytes'   : 1024 * 1024 * 5,                              # 5MB
            'backupCount': 5,
            'formatter'  : 'standard',
        },
        'computations': {
            'class': 'logging.handlers.RotatingFileHandler',
            'filename': os.path.join(BASE_DIR,'logs','computations.log'),
            'maxBytes': 1024 * 1024 * 5, # 5MB
            'backupCount': 5,
            'formatter': 'standard',
        },
    
    },
    'loggers': {
        'general': {
            'handlers': ['general'],
            'level': 'DEBUG',
        },
        'updates': {
            'handlers': ['updates'],
            'level': 'DEBUG',
        },
        'accesses': {
            'handlers': ['accesses'],
            'level': 'DEBUG',
        },
        'computations': {
            'handlers': ['computations'],
            'level': 'DEBUG',
        },
    },
}

STATIC_URL = '/static/'
STATIC_ROOT = os.path.join(BASE_DIR, 'staticfiles')
STATICFILES_DIRS = [
    os.path.join(BASE_DIR, 'static'),
]
