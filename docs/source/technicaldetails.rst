Technical Details
=================

Connection Recycling
---------------------
Building the connection to the remote host is a computationally expensive operation, which is why lwreg caches the connection to the database as default behaviour.
Whenever you work from lwreg and are just using a single connection, making use of the caching is the most efficient approach.
However, in a scenario where you are working with multiple connections to the same database (for example if you have a lot of jobs that all use the database running), it is usually better to disable the caching in order to avoid exceeding the databases's connection limit.
When calling :code:`configure_from_database`, there is a function argument called :code:`cacheConnection` which allows you to set the caching behaviour to :code:`False`.
This will disable the caching and create a new connection to the database every time you call :code:`configure_from_database`.

For the rest of the lwreg functions, the caching behaviour is set through the passed :code:`config` setting of :code:`cacheConnection`.