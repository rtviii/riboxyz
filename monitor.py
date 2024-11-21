from datetime import datetime
import psutil
import docker
import json
from pathlib import Path
import logging

class ResourceMonitor:
    def __init__(self):
        self.docker_client = docker.from_env()
        self.metrics_dir = Path("./metrics")
        self.metrics_dir.mkdir(parents=True, exist_ok=True)

    def get_host_metrics(self):
        """Get host system metrics"""
        return {
            "memory": {
                "total_gb": psutil.virtual_memory().total / (1024**3),
                "available_gb": psutil.virtual_memory().available / (1024**3),
                "percent_used": psutil.virtual_memory().percent
            },
            "cpu": {
                "percent": psutil.cpu_percent(interval=1),
                "load_avg": psutil.getloadavg(),
                "core_count": psutil.cpu_count()
            },
            "disk": {
                "root": {
                    "total_gb"    : psutil.disk_usage('/home/rtviii/dev').total   / (1024**3),
                    "free_gb"     : psutil.disk_usage('/home/rtviii/dev').free    / (1024**3),
                    "percent_used": psutil.disk_usage('/home/rtviii/dev').percent
                },
                "data": {
                    "total_gb"    : psutil.disk_usage('/home/rtviii/dev/ribxz_neo4j_volume/data').total   / (1024**3),
                    "free_gb"     : psutil.disk_usage('/home/rtviii/dev/ribxz_neo4j_volume/data').free    / (1024**3),
                    "percent_used": psutil.disk_usage('/home/rtviii/dev/ribxz_neo4j_volume/data').percent
                }
            }
        }

    def get_container_metrics(self):
        """Get container-specific metrics"""
        metrics = {}
        for container in self.docker_client.containers.list():
            if container.name in ['django', 'neo']:
                stats = container.stats(stream=False)
                
                # Get memory limits if set
                mem_limit = stats['memory_stats'].get('limit', 0)
                mem_usage = stats['memory_stats'].get('usage', 0)
                
                # Calculate CPU usage
                cpu_delta = stats['cpu_stats']['cpu_usage']['total_usage'] - \
                           stats['precpu_stats']['cpu_usage']['total_usage']
                system_delta = stats['cpu_stats']['system_cpu_usage'] - \
                             stats['precpu_stats']['system_cpu_usage']
                cpu_percent = 0.0
                if system_delta > 0:
                    cpu_percent = (cpu_delta / system_delta) * len(stats['cpu_stats']['cpu_usage'].get('percpu_usage', [1])) * 100

                metrics[container.name] = {
                    "status": container.status,
                    "memory": {
                        "usage_gb": mem_usage / (1024**3),
                        "limit_gb": mem_limit / (1024**3) if mem_limit else "unlimited",
                        "percent": (mem_usage / mem_limit * 100) if mem_limit else None
                    },
                    "cpu": {
                        "percent": round(cpu_percent, 2),
                        "throttling_data": stats['cpu_stats'].get('throttling_data', {})
                    },
                    "restarts": container.attrs['RestartCount']
                }
        return metrics

    def check_thresholds(self, metrics):
        """Check for concerning resource usage"""
        warnings = []
        
        # Host memory warnings
        if metrics['host']['memory']['percent_used'] > 85:
            warnings.append(f"High host memory usage: {metrics['host']['memory']['percent_used']}%")
        
        # Host disk warnings
        if metrics['host']['disk']['root']['percent_used'] > 85:
            warnings.append(f"High root disk usage: {metrics['host']['disk']['root']['percent_used']}%")
        if metrics['host']['disk']['data']['percent_used'] > 85:
            warnings.append(f"High data disk usage: {metrics['host']['disk']['data']['percent_used']}%")

        # Container warnings
        for container, stats in metrics['containers'].items():
            if stats['status'] != 'running':
                warnings.append(f"Container {container} is {stats['status']}")
            if stats['memory'].get('percent') and stats['memory']['percent'] > 90:
                warnings.append(f"High memory usage in {container}: {stats['memory']['percent']}%")
            if stats['restarts'] > 0:
                warnings.append(f"Container {container} has restarted {stats['restarts']} times")

        return warnings

    def collect_metrics(self):
        """Collect all metrics"""
        metrics = {
            "timestamp": datetime.now().isoformat(),
            "host": self.get_host_metrics(),
            "containers": self.get_container_metrics()
        }
        
        warnings = self.check_thresholds(metrics)
        if warnings:
            metrics["warnings"] = warnings

        # Save metrics
        metrics_file = self.metrics_dir / f"metrics_{datetime.now().strftime('%Y%m%d')}.json"
        with open(metrics_file, 'a') as f:
            json.dump(metrics, f)
            f.write('\n')

        return metrics

if __name__ == "__main__":
    monitor = ResourceMonitor()
    metrics = monitor.collect_metrics()
    print(json.dumps(metrics, indent=2))