Examples
========

Bulk registration and duplicates
--------------------------------

The following example session demonstrates how to register multiple compounds at once and how to change what happens if duplicates are present.

We start by bulk registering all of the molecules in an SDF and looking at the molregnos returned::
    
  >>> mrns = lwreg.bulk_register(sdfile='./mols1.sdf')
  >>> mrns
  (1, 2, 3)

Now we bulk register a set of molecules from a different file that has two duplicates (molecules already present in the database) and can see that the result shows that the duplicates were the second an third molecules in the file::

  >>> mrns = lwreg.bulk_register(sdfile='./mols2.sdf')
  >>> mrns
  (4, <RegistrationFailureReasons.DUPLICATE: 1>, <RegistrationFailureReasons.DUPLICATE: 1>)


If we use the ``fail_on_duplicate`` parameter, we can change the behavior of the bulk registration to return the molregnos of the duplicates instead of the ``RegistrationFailureReasons.DUPLICATE`` marker::

  >>> mrns = lwreg.bulk_register(sdfile='./mols2.sdf', fail_on_duplicate=False)
  >>> mrns
  (4, 3, 1)


