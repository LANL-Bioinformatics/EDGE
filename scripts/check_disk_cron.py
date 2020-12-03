from crontab import CronTab
import os
cron = CronTab(user='l197233')
cron.remove_all()
print(os.getcwd())
job = cron.new(command='/Users/l197233/anaconda3/envs/edge-monitor/bin/python /Devel/edge/scripts/check_disk_usage.py -v /Users -t 0.8 >> /Devel/edge/scripts/disk_usage.log 2>&1')
job.env['PATH'] = '/Users/l197233/anaconda3/envs/edge-monitor/bin/python'
job.minute.every(1)

cron.write()