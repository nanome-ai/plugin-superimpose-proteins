import os
from cookbook.interface import PluginInstanceRedisInterface
# Set up redis credentials
redis_host = 'redis'
redis_port = 6379
redis_password = ''
redis_channel = os.environ.get('REDIS_CHANNEL')

# When your PluginService is running, you can get your channel value from the Logs, or from the query parameter in your open browser.
# Update this value to match that, so that your commands will run against your live workspace.
redis_channel = os.environ.get("REDIS_CHANNEL")

# Increase the recursion limit in order to properly serialize Complexes
# recursion_limit = 10000
# sys.setrecursionlimit(recursion_limit)

plugin_instance = PluginInstanceRedisInterface(redis_host, redis_port, redis_password, redis_channel)
comps = plugin_instance.request_complex_list()