#! /bin/bash
uvicorn main:app --host 192.168.2.233 --port 5050 --log-config log_config.ini
