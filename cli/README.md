oclif-hello-world
=================

oclif example Hello World CLI

[![oclif](https://img.shields.io/badge/cli-oclif-brightgreen.svg)](https://oclif.io)
[![Version](https://img.shields.io/npm/v/oclif-hello-world.svg)](https://npmjs.org/package/oclif-hello-world)
[![CircleCI](https://circleci.com/gh/oclif/hello-world/tree/main.svg?style=shield)](https://circleci.com/gh/oclif/hello-world/tree/main)
[![Downloads/week](https://img.shields.io/npm/dw/oclif-hello-world.svg)](https://npmjs.org/package/oclif-hello-world)
[![License](https://img.shields.io/npm/l/oclif-hello-world.svg)](https://github.com/oclif/hello-world/blob/main/package.json)

<!-- toc -->
* [Usage](#usage)
* [Commands](#commands)
<!-- tocstop -->
# Usage
<!-- usage -->
```sh-session
$ npm install -g ribxzcli
$ ribxzcli COMMAND
running command...
$ ribxzcli (--version)
ribxzcli/0.3.0 linux-x64 node-v18.14.0
$ ribxzcli --help [COMMAND]
USAGE
  $ ribxzcli COMMAND
...
```
<!-- usagestop -->
# Commands
<!-- commands -->
* [`ribxzcli db`](#ribxzcli-db)
* [`ribxzcli db status`](#ribxzcli-db-status)
* [`ribxzcli help [COMMAND]`](#ribxzcli-help-command)
* [`ribxzcli struct RCSB_ID`](#ribxzcli-struct-rcsb_id)
* [`ribxzcli struct obtain RCSB_ID`](#ribxzcli-struct-obtain-rcsb_id)
* [`ribxzcli struct show RCSB_ID`](#ribxzcli-struct-show-rcsb_id)

## `ribxzcli db`

```
USAGE
  $ ribxzcli db [-r <value>] [-a <value>] [-p <value>] [-u <value>] [-d <value>] [-e <value>]

FLAGS
  -a, --NEO4J_URI=<value>
  -d, --NEO4J_CURRENTDB=<value>
  -e, --env=<value>              Environment variable
  -p, --NEO4J_PASSWORD=<value>
  -r, --RIBETL_DATA=<value>
  -u, --NEO4J_USER=<value>
```

_See code: [dist/commands/db/index.js](https://github.com/rtviii/hello-world/blob/v0.3.0/dist/commands/db/index.js)_

## `ribxzcli db status`

Query structure in the database

```
USAGE
  $ ribxzcli db status [-r <value>] [-a <value>] [-p <value>] [-u <value>] [-d <value>] [--PYTHONBIN <value>]
    [--COMMIT_STRUCTURE_SH <value>] [--EXTRACT_BSITES_PY <value>] [--RENDER_THUMBNAIL_PY <value>] [--SPLIT_RENAME_PY
    <value>] [-e <value>] [-A]

FLAGS
  -A, --updateall
  -a, --NEO4J_URI=<value>
  -d, --NEO4J_CURRENTDB=<value>
  -e, --env=<value>              Environment variable
  -p, --NEO4J_PASSWORD=<value>
  -r, --RIBETL_DATA=<value>
  -u, --NEO4J_USER=<value>
  --COMMIT_STRUCTURE_SH=<value>
  --EXTRACT_BSITES_PY=<value>
  --PYTHONBIN=<value>
  --RENDER_THUMBNAIL_PY=<value>
  --SPLIT_RENAME_PY=<value>

DESCRIPTION
  Query structure in the database
```

## `ribxzcli help [COMMAND]`

Display help for ribxzcli.

```
USAGE
  $ ribxzcli help [COMMAND] [-n]

ARGUMENTS
  COMMAND  Command to show help for.

FLAGS
  -n, --nested-commands  Include all nested commands in the output.

DESCRIPTION
  Display help for ribxzcli.
```

_See code: [@oclif/plugin-help](https://github.com/oclif/plugin-help/blob/v5.1.22/src/commands/help.ts)_

## `ribxzcli struct RCSB_ID`

```
USAGE
  $ ribxzcli struct [RCSB_ID] [-r <value>] [-a <value>] [-p <value>] [-u <value>] [-d <value>] [-e <value>]

FLAGS
  -a, --NEO4J_URI=<value>
  -d, --NEO4J_CURRENTDB=<value>
  -e, --env=<value>              Environment variable
  -p, --NEO4J_PASSWORD=<value>
  -r, --RIBETL_DATA=<value>
  -u, --NEO4J_USER=<value>
```

_See code: [dist/commands/struct/index.js](https://github.com/rtviii/hello-world/blob/v0.3.0/dist/commands/struct/index.js)_

## `ribxzcli struct obtain RCSB_ID`

Query structure in the database

```
USAGE
  $ ribxzcli struct obtain [RCSB_ID] [-r <value>] [-a <value>] [-p <value>] [-u <value>] [-d <value>] [--PYTHONBIN
    <value>] [--COMMIT_STRUCTURE_SH <value>] [--EXTRACT_BSITES_PY <value>] [--RENDER_THUMBNAIL_PY <value>]
    [--SPLIT_RENAME_PY <value>] [-e <value>] [-R] [-C]

FLAGS
  -C, --commit
  -R, --repair
  -a, --NEO4J_URI=<value>
  -d, --NEO4J_CURRENTDB=<value>
  -e, --env=<value>              Environment variable
  -p, --NEO4J_PASSWORD=<value>
  -r, --RIBETL_DATA=<value>
  -u, --NEO4J_USER=<value>
  --COMMIT_STRUCTURE_SH=<value>
  --EXTRACT_BSITES_PY=<value>
  --PYTHONBIN=<value>
  --RENDER_THUMBNAIL_PY=<value>
  --SPLIT_RENAME_PY=<value>

DESCRIPTION
  Query structure in the database
```

## `ribxzcli struct show RCSB_ID`

Query structure in the database

```
USAGE
  $ ribxzcli struct show [RCSB_ID] [-r <value>] [-a <value>] [-p <value>] [-u <value>] [-d <value>] [--PYTHONBIN
    <value>] [--COMMIT_STRUCTURE_SH <value>] [--EXTRACT_BSITES_PY <value>] [--RENDER_THUMBNAIL_PY <value>]
    [--SPLIT_RENAME_PY <value>] [-e <value>] [--files] [--db] [--dryrun] [-f]

FLAGS
  -a, --NEO4J_URI=<value>
  -d, --NEO4J_CURRENTDB=<value>
  -e, --env=<value>              Environment variable
  -f, --force
  -p, --NEO4J_PASSWORD=<value>
  -r, --RIBETL_DATA=<value>
  -u, --NEO4J_USER=<value>
  --COMMIT_STRUCTURE_SH=<value>
  --EXTRACT_BSITES_PY=<value>
  --PYTHONBIN=<value>
  --RENDER_THUMBNAIL_PY=<value>
  --SPLIT_RENAME_PY=<value>
  --db
  --dryrun
  --files

DESCRIPTION
  Query structure in the database
```
<!-- commandsstop -->
