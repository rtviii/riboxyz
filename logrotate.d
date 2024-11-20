# /etc/logrotate.d/ribxz
/var/log/ribxz/*.log {
    rotate 5
    size 50M
    compress
    delaycompress
    missingok
    notifempty
    create 0644 rtviii docker
    postrotate
        systemctl kill -s HUP ribxz-update.service || true
    endscript
}