Lexer
=====

Get the Code
------------
Get the pygment source code from bitbucket repo with :code:`hg clone http://bitbucket.org/birkenfeld/pygments-main pygments`

Register Your Lexer
-------------------
To make Pygments aware of your new lexer, you have to perform the following steps:

1. First, change to the current directory containing the Pygments source code:

    .. code::

        $ cd .../pygments

   Select a matching module under pygments/lexers, or create a new module for your lexer class.
   Copy the BELLexer.py in that module


2. The lexer can be made publicly known by rebuilding the lexer mapping:

    .. code::

        $ make mapfiles

Test
----
To test the new lexer, store an example file with the proper extension in tests/examplefiles.
For example, to test the BELLexer, add a tests/examplefiles/Example.bel containing a sample bel code.

Now you can use pygmentize to render your example to HTML:

.. code::

    $ ./pygmentize -O full -f html -o /tmp/example.html tests/examplefiles/Example.bel

To view the result, open ./example.html in your browser.
