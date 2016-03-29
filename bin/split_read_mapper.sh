#!/bin/bash

BIN="$( dirname $( which $0 ) )"
WT_HOME="$( dirname $BIN )"
TERM_WIDTH="$(stty size | sed 's/.* //;')"

echo USING WT_HOME=$WT_HOME 1>&2

exec java -Dcom.lifetechnologies.solid.wt.home="$WT_HOME" -Dterm.width="$TERM_WIDTH" -Dcommand="$0" -cp "$WT_HOME/pkg/WholeTranscriptome.jar:$WT_HOME/lib/*" com.lifetechnologies.solid.wt.cmd.CmdRunner com.lifetechnologies.solid.wt.cmd.MapperCmd $@
