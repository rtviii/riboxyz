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
* [`ribxz hello PERSON`](#ribxz-hello-person)
* [`ribxz hello world`](#ribxz-hello-world)
* [`ribxz help [COMMAND]`](#ribxz-help-command)
* [`ribxz plugins`](#ribxz-plugins)
* [`ribxz plugins:install PLUGIN...`](#ribxz-pluginsinstall-plugin)
* [`ribxz plugins:inspect PLUGIN...`](#ribxz-pluginsinspect-plugin)
* [`ribxz plugins:install PLUGIN...`](#ribxz-pluginsinstall-plugin-1)
* [`ribxz plugins:link PLUGIN`](#ribxz-pluginslink-plugin)
* [`ribxz plugins:uninstall PLUGIN...`](#ribxz-pluginsuninstall-plugin)
* [`ribxz plugins:uninstall PLUGIN...`](#ribxz-pluginsuninstall-plugin-1)
* [`ribxz plugins:uninstall PLUGIN...`](#ribxz-pluginsuninstall-plugin-2)
* [`ribxz plugins update`](#ribxz-plugins-update)

## `ribxz hello PERSON`

Say hello

```
USAGE
  $ ribxz hello [PERSON] -f <value>

ARGUMENTS
  PERSON  Person to say hello to

FLAGS
  -f, --from=<value>  (required) Who is saying hello

DESCRIPTION
  Say hello

EXAMPLES
  $ oex hello friend --from oclif
  hello friend from oclif! (./src/commands/hello/index.ts)
```

_See code: [dist/commands/hello/index.ts](https://github.com/rtviii/riboxyz/blob/v0.0.0/dist/commands/hello/index.ts)_

## `ribxz hello world`

Say hello world

```
USAGE
  $ ribxz hello world

DESCRIPTION
  Say hello world

EXAMPLES
  $ ribxz hello world
  hello world! (./src/commands/hello/world.ts)
```

## `ribxz help [COMMAND]`

Display help for ribxz.

```
USAGE
  $ ribxz help [COMMAND] [-n]

ARGUMENTS
  COMMAND  Command to show help for.

FLAGS
  -n, --nested-commands  Include all nested commands in the output.

DESCRIPTION
  Display help for ribxz.
```

_See code: [@oclif/plugin-help](https://github.com/oclif/plugin-help/blob/v5.1.19/src/commands/help.ts)_

## `ribxz plugins`

List installed plugins.

```
USAGE
  $ ribxz plugins [--core]

FLAGS
  --core  Show core plugins.

DESCRIPTION
  List installed plugins.

EXAMPLES
  $ ribxz plugins
```

_See code: [@oclif/plugin-plugins](https://github.com/oclif/plugin-plugins/blob/v2.1.7/src/commands/plugins/index.ts)_

## `ribxz plugins:install PLUGIN...`

Installs a plugin into the CLI.

```
USAGE
  $ ribxz plugins:install PLUGIN...

ARGUMENTS
  PLUGIN  Plugin to install.

FLAGS
  -f, --force    Run yarn install with force flag.
  -h, --help     Show CLI help.
  -v, --verbose

DESCRIPTION
  Installs a plugin into the CLI.
  Can be installed from npm or a git url.

  Installation of a user-installed plugin will override a core plugin.

  e.g. If you have a core plugin that has a 'hello' command, installing a user-installed plugin with a 'hello' command
  will override the core plugin implementation. This is useful if a user needs to update core plugin functionality in
  the CLI without the need to patch and update the whole CLI.


ALIASES
  $ ribxz plugins add

EXAMPLES
  $ ribxz plugins:install myplugin 

  $ ribxz plugins:install https://github.com/someuser/someplugin

  $ ribxz plugins:install someuser/someplugin
```

## `ribxz plugins:inspect PLUGIN...`

Displays installation properties of a plugin.

```
USAGE
  $ ribxz plugins:inspect PLUGIN...

ARGUMENTS
  PLUGIN  [default: .] Plugin to inspect.

FLAGS
  -h, --help     Show CLI help.
  -v, --verbose

DESCRIPTION
  Displays installation properties of a plugin.

EXAMPLES
  $ ribxz plugins:inspect myplugin
```

## `ribxz plugins:install PLUGIN...`

Installs a plugin into the CLI.

```
USAGE
  $ ribxz plugins:install PLUGIN...

ARGUMENTS
  PLUGIN  Plugin to install.

FLAGS
  -f, --force    Run yarn install with force flag.
  -h, --help     Show CLI help.
  -v, --verbose

DESCRIPTION
  Installs a plugin into the CLI.
  Can be installed from npm or a git url.

  Installation of a user-installed plugin will override a core plugin.

  e.g. If you have a core plugin that has a 'hello' command, installing a user-installed plugin with a 'hello' command
  will override the core plugin implementation. This is useful if a user needs to update core plugin functionality in
  the CLI without the need to patch and update the whole CLI.


ALIASES
  $ ribxz plugins add

EXAMPLES
  $ ribxz plugins:install myplugin 

  $ ribxz plugins:install https://github.com/someuser/someplugin

  $ ribxz plugins:install someuser/someplugin
```

## `ribxz plugins:link PLUGIN`

Links a plugin into the CLI for development.

```
USAGE
  $ ribxz plugins:link PLUGIN

ARGUMENTS
  PATH  [default: .] path to plugin

FLAGS
  -h, --help     Show CLI help.
  -v, --verbose

DESCRIPTION
  Links a plugin into the CLI for development.
  Installation of a linked plugin will override a user-installed or core plugin.

  e.g. If you have a user-installed or core plugin that has a 'hello' command, installing a linked plugin with a 'hello'
  command will override the user-installed or core plugin implementation. This is useful for development work.


EXAMPLES
  $ ribxz plugins:link myplugin
```

## `ribxz plugins:uninstall PLUGIN...`

Removes a plugin from the CLI.

```
USAGE
  $ ribxz plugins:uninstall PLUGIN...

ARGUMENTS
  PLUGIN  plugin to uninstall

FLAGS
  -h, --help     Show CLI help.
  -v, --verbose

DESCRIPTION
  Removes a plugin from the CLI.

ALIASES
  $ ribxz plugins unlink
  $ ribxz plugins remove
```

## `ribxz plugins:uninstall PLUGIN...`

Removes a plugin from the CLI.

```
USAGE
  $ ribxz plugins:uninstall PLUGIN...

ARGUMENTS
  PLUGIN  plugin to uninstall

FLAGS
  -h, --help     Show CLI help.
  -v, --verbose

DESCRIPTION
  Removes a plugin from the CLI.

ALIASES
  $ ribxz plugins unlink
  $ ribxz plugins remove
```

## `ribxz plugins:uninstall PLUGIN...`

Removes a plugin from the CLI.

```
USAGE
  $ ribxz plugins:uninstall PLUGIN...

ARGUMENTS
  PLUGIN  plugin to uninstall

FLAGS
  -h, --help     Show CLI help.
  -v, --verbose

DESCRIPTION
  Removes a plugin from the CLI.

ALIASES
  $ ribxz plugins unlink
  $ ribxz plugins remove
```

## `ribxz plugins update`

Update installed plugins.

```
USAGE
  $ ribxz plugins update [-h] [-v]

FLAGS
  -h, --help     Show CLI help.
  -v, --verbose

DESCRIPTION
  Update installed plugins.
```
<!-- commandsstop -->
