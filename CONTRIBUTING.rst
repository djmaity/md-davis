Contributing
============

Contributions are welcome, and they are greatly appreciated! Every little bit
helps, and credit will always be given.

Types of Contributions
----------------------

Report Bugs
^^^^^^^^^^^

When `reporting a bug <https://github.com/djmaity/md-davis/issues>`_ please include:

* Your operating system name and version.
* Any details about your local setup that might be helpful in troubleshooting.
* Detailed steps to reproduce the bug.

Fix Bugs
^^^^^^^^

Look through the GitHub issues for bugs. Anything tagged with "bug" and "help
wanted" is open to whoever wants to implement it.

Implement Features
^^^^^^^^^^^^^^^^^^

Look through the GitHub issues for features. Anything tagged with "enhancement"
and "help wanted" is open to whoever wants to implement it.

Write Documentation
^^^^^^^^^^^^^^^^^^^

MD DaVis could always use more documentation, whether as part of the
official MD DaVis docs, in docstrings, or even on the web in blog posts,
articles, and such.

Feature Requests and Feedback
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The best way to send feedback is to file an `issue <https://github.com/djmaity/md-davis/issues>`_.

If you are proposing a feature:

* Explain in detail how it would work.
* Keep the scope as narrow as possible, to make it easier to implement.
* Remember that this is a volunteer-driven project, and that contributions
  are welcome :)


Development
-----------
..
    To run all the tests run::

        tox

    Note, to combine the coverage data from all the tox environments run:

    .. list-table::
        :widths: 10 90
        :stub-columns: 1

        - - Windows
          - ::

                set PYTEST_ADDOPTS=--cov-append
                tox

        - - Other
          - ::

                PYTEST_ADDOPTS=--cov-append tox


    Development Environment
    ^^^^^^^^^^^^^^^^^^^^^^^

Create the development environment and install the dependencies using the
`dev_environment.yml <https://github.com/djmaity/md-davis/blob/master/dev_environment.yml>`_ file.
This installs packages for linting, packaging and building documentation in addition to the core dependencies.

.. code-block:: bash

    conda env create -f dev_environment.yml -n md_davis_dev
    conda activate md_davis_dev

Install md-davis in editable mode:

.. code-block:: bash

    pip install -e md-davis

..
    To set up `md_davis` for local development:

    1. Fork `md_davis <https://github.com/djmaity/md_davis>`_ (look for the "Fork" button).

    2. Clone your fork locally::

        git clone git@github.com:YOURGITHUBNAME/md_davis.git

    3. Install your local copy into a virtualenv. Assuming you have virtualenvwrapper installed, this is how you set up your fork for local development::

        $ mkvirtualenv md_davis
        $ cd md_davis/
        $ python setup.py develop

    4. Create a branch for local development::

        git checkout -b name-of-your-bugfix-or-feature

       Now you can make your changes locally.

    5. When you're done making changes, check that your changes pass flake8, doc builder and the
       tests, including testing other Python versions with tox::

        $ flake8 md_davis tests
        $ python setup.py test or py.test
        $ tox

       To get flake8 and tox, just pip install them into your virtualenv.

    6. Commit your changes and push your branch to GitHub::

        git add .
        git commit -m "Your detailed description of your changes."
        git push origin name-of-your-bugfix-or-feature

    7. Submit a pull request through the GitHub website.


    Pull Request Guidelines
    ^^^^^^^^^^^^^^^^^^^^^^^

    If you need some code review or feedback while you're developing the code
    just make the pull request. Before you submit a pull request, check that
    it meets these guidelines:

    1. The pull request should include passing tests (run ``tox``).

    2. If the pull request adds functionality, the docs should be updated. Put
       your new functionality into a function with a docstring, and add the
       feature to the list in README.rst.

    3. The pull request should work for Python 3.7, 3.8 and 3.9, and for PyPy. Check
       https://travis-ci.org/djmaity/md_davis/pull_requests
       and make sure that the tests pass for all supported Python versions.

    4. Add a note to HISTORY.rst about the changes.

    5. Add yourself to authors in README.md.

    Tips
    ^^^^

    To run a subset of tests::

        tox -e envname -- pytest -k test_myfeature

    To run all the test environments in *parallel*::

        tox -p auto

    Deploying
    ^^^^^^^^^

    A reminder for the maintainers on how to deploy.
    Make sure all your changes are committed (including an entry in HISTORY.rst).
    Then run::

    $ bumpversion patch     # possible: major / minor / patch
    $ git push
    $ git push --tags

    Travis will then deploy to PyPI if tests pass.

