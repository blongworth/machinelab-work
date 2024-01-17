#!/bin/bash

# render and publish LECS data from website
# use data from the last week
# publish to https://blongworth.github.io/machinelab-work/

quarto render read_parse_web.qmd -P start_date:$(date -v-7d "+%Y-%m-%d")
quarto publish read_parse_web.qmd --no-render --no-prompt --no-browser

