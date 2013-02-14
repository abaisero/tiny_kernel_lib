Tutorial to Markdown {#tutorial}
====================

Level 1 Header
==============

Level 2 Header
--------------

Block quotes
------------

> This is a block quote
> spanning multiple lines

Lists
-----

- item 1
- item 2
  + nested item 1
  + nested item 2
    * again..
    * and again!
- item 3

Code Blocks
-----------

    int a=5;

Horizontal Rulers
-----------------

---

___

***

Emphasis
--------

*single asterisks*

_single underscores_

**double asterisks**

__double underscores__

Inline Code Block
-----------------

Like, for example, `int a=5;`

Links
-----

[Link to Google, w/o overtext](http://www.google.it)

[Link to Google, w/ overtext](http://www.google.it "overtext")

[Link to RefKernel class](@ref RefKernel)

[link_to_google]: http://www.google.it "overtext"

[Non inline link][link_to_google]

[link_to_google][]

Images
------

Let's skip that for now..

Automatic Link
--------------

<http://www.google.it>

Tables
------

A     | Simple
------|--------
Table | Example
With  | Multiple
Rows  | (No empty cells allowed)

Fenced Code Blocks
------------------

~~~~~~~~~~{.cpp}
// kinda nicer than the other code blocks
int a=5;
int b=5;
int c=a+b;
~~~~~~~~~~

