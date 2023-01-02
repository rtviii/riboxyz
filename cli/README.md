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
$ npm install -g ribxz
$ ribxz COMMAND
running command...
$ ribxz (--version)
ribxz/0.0.0 linux-x64 node-v14.17.6
$ ribxz --help [COMMAND]
USAGE
  $ ribxz COMMAND
...
```
<!-- usagestop -->
# Commands
<!-- commands -->
* [`ribxz db`](#ribxz-db)
* [`ribxz struct RCSB_ID`](#ribxz-struct-rcsb_id)
* [`ribxz struct show RCSB_ID`](#ribxz-struct-show-rcsb_id)

## `ribxz db`

describe the command here

```
USAGE
  $ ribxz db [-r <value>] [-a <value>] [-p <value>] [-u <value>] [-d <value>] [--PYTHONBIN <value>]
    [--COMMIT_STRUCTURE_SH <value>] [--EXTRACT_BSITES_PY <value>] [--RENDER_THUMBNAIL_PY <value>] [--SPLIT_RENAME_PY
    <value>] [-e <value>] [-f]

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

DESCRIPTION
  describe the command here
```

_See code: [dist/commands/db/index.js](https://github.com/rtviii/riboxyz/blob/v0.0.0/dist/commands/db/index.js)_

## `ribxz struct RCSB_ID`

```
USAGE
  $ ribxz struct [RCSB_ID] [-r <value>] [-a <value>] [-p <value>] [-u <value>] [-d <value>] [-e <value>]

FLAGS
  -a, --NEO4J_URI=<value>
  -d, --NEO4J_CURRENTDB=<value>
  -e, --env=<value>              Environment variable
  -p, --NEO4J_PASSWORD=<value>
  -r, --RIBETL_DATA=<value>
  -u, --NEO4J_USER=<value>
```

_See code: [dist/commands/struct/index.js](https://github.com/rtviii/riboxyz/blob/v0.0.0/dist/commands/struct/index.js)_

## `ribxz struct show RCSB_ID`

Query structure in the database

```
USAGE
  $ ribxz struct show [RCSB_ID] [-r <value>] [-a <value>] [-p <value>] [-u <value>] [-d <value>] [--PYTHONBIN
    <value>] [--COMMIT_STRUCTURE_SH <value>] [--EXTRACT_BSITES_PY <value>] [--RENDER_THUMBNAIL_PY <value>]
    [--SPLIT_RENAME_PY <value>] [-e <value>] [--files] [--db] [-R] [-C] [--dryrun] [-f]

FLAGS
  -C, --commit
  -R, --repair
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
