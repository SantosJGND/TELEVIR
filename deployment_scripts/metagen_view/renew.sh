#!/bin/bash


read -p "Press enter to continue"

rm db.sqlite3
python manage.py migrate
