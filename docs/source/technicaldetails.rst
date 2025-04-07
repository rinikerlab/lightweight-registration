Technical Details
=================

Connection Recycling
---------------------
Building the connection to the remote host is a computationally expensive operation, which is why lwreg caches the connection to the database as default behaviour.
Whenever you work from lwreg connecting to the same database with just a single connection, making use of the caching is the most efficient.
However, in a scenario where you are working with multiple connections to the same database, it is advised to disable the caching.
When calling :code:`configure_from_database`, there is a function argument called :code:`cache_connection` which allows you to set the caching behaviour to :code:`False`.
This will disable the caching and create a new connection to the database every time you call :code:`configure_from_database`.

For the rest of the lwreg functions, the caching behaviour is set through the passed :code:`config` setting of :code:`cacheConnection`.