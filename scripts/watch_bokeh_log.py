from log_watcher import LogWatcher
from send_notifications import send_sms_message, send_email_message
import os
import shutil
import time

from datetime import datetime
def callback(filename, lines):
    print(filename, lines)
    for line in lines:
        if 'NoneType' in str(line):
            print('Bokeh server restarting')
            shutil.copyfile('/panfs/biopan01/edge_prod/edge_ui/bokeh.log',
                            '/panfs/biopan01/edge_prod/edge_ui/logs/bokeh.log-{0}'.format(datetime.now().isoformat()))
            send_sms_message('Bokeh server restarting: {0}'.format(str(line)))
            send_email_message('Bokeh server restarting: {0}'.format(str(line)))
            os.system('/panfs/biopan01/edge_prod/scripts/start_bokeh.sh')
            break
lw = LogWatcher('/panfs/biopan01/edge_prod/edge_ui/', callback, tail_lines=10)
idx = 0
while 1:
    lw.loop(blocking=False)
    time.sleep(0.1)
    idx += 1
    if idx % 1000 == 0:
        print('Checking server...{}'.format(datetime.now().isoformat()))