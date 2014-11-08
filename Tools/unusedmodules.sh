#!/bin/bash
comm -2 -3 <(grep -rl "use $1" *|sort) <(grep -rl "\<$2" *|sort)
