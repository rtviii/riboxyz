import os


# SECURITY WARNING: keep the secret key used in production secret!
SECRET_KEY = 'ju=n4om3z00jd1+y2(ufn)g^@w-dj*&-45&4yd1_aiun50b6by'

BASE_DIR     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RIBETL_DATA  = os.path.join(BASE_DIR, "ribetldata"                           ) # this should be either docker-mounted or populated through the utils module
CYPHER_EXEC  = os.path.join(BASE_DIR, "mod_db","cypher_ops","cypher_exec")     # this should be either docker-mounted or populated through the utils module

neo4j_vars = ["NEO4J_URI", "NEO4J_USER", "NEO4J_PASSWORD", "NEO4J_CURRENTDB"]
print(" ---------------------------------------------- App has been reset. -----------------------------------------------") 
for var in neo4j_vars:
    # print("Environment variable {}:\t{}".format( var, os.getenv(var) ))
    if var not in os.environ:
        print("Environment variable {} not set".format(var))
        exit(1)

NEO4J_URI       = os.getenv("NEO4J_URI")
NEO4J_PASSWORD  = os.getenv("NEO4J_PASSWORD")
NEO4J_USER      = os.getenv("NEO4J_USER")
NEO4J_CURRENTDB = os.getenv("NEO4J_CURRENTDB")

DEBUG         = True
ALLOWED_HOSTS = ["*"]

INSTALLED_APPS = [
    'django.contrib.admin'       ,
    'django.contrib.auth'        ,
    'django.contrib.contenttypes',
    'django.contrib.sessions'    ,
    'django.contrib.messages'    ,
    'django.contrib.staticfiles' ,
    'rest_framework',
    'mod_comp',
    'mod_db',
    'mod_utils',
]

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
        'BACKEND' : 'django.template.backends.django.DjangoTemplates',
        'DIRS'    : [],
        'APP_DIRS': True,
        'OPTIONS' : {
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


# Database
# https://docs.djangoproject.com/en/1.11/ref/settings/#databases

DATABASES = {
    'default': {
        'ENGINE': 'django.db.backends.sqlite3',
        'NAME'  : os.path.join(BASE_DIR, 'db.sqlite3'),
    }
}


# Password validation
# https://docs.djangoproject.com/en/1.11/ref/settings/#auth-password-validators

AUTH_PASSWORD_VALIDATORS = [
    {'NAME': 'django.contrib.auth.password_validation.UserAttributeSimilarityValidator',},
    {'NAME': 'django.contrib.auth.password_validation.MinimumLengthValidator',},
    {'NAME': 'django.contrib.auth.password_validation.CommonPasswordValidator',},
    {'NAME': 'django.contrib.auth.password_validation.NumericPasswordValidator',},
]


# Internationalization
# https://docs.djangoproject.com/en/1.11/topics/i18n/

LANGUAGE_CODE    = 'en-us'
TIME_ZONE        = 'UTC'
USE_I18N         = True
USE_L10N         = True
USE_TZ           = True
STATIC_URL       = '/static/'
STATIC_ROOT      = '/static/'
# # STATIC_ROOT      = os.path.join(BASE_DIR, 'rbxz_bend/static')
STATICFILES_DIRS = (
    '/static/',
    # os.path.join(BASE_DIR, 'rbxz_bend/static'),
)

