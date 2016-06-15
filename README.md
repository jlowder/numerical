# numerical

A library of mumerical methods for Common Lisp. Currently a work in progress.

## Overview

The library consists of mostly higher-order functions, i.e. functions
that take a function as a parameter and return a function as the
result. For example,

~~~lisp
(ql:quickload :numerical)
(use-package '(:numerical :numerical.roots))

(defun f (x)
  (sin x))

(bisect #'f)
~~~

This will create a closure that uses the "bisect" method of root
finding on the function "f". Each call to the closure will return 3
values: the numerical result, a boolean indicating whether or not it
converged on the solution to within tolerance, and the number of
iteration attempts performed.

~~~lisp
(funcall (bisect #'f) 1 5)

3.1415926814079285d0
T
25

(funcall (bisect #'f) 1 3)
2.0d0
NIL
0
~~~

The default tolerance is 8 decimal places (5d-8) and the default iteration limit is 50. These can be changed with "tolerance" and "limit" keyword parameters:

~~~lisp
(funcall (bisect #'f :tolerance 5d-12 :limit 100) 1 5)

3.1415926535919425d0
T
38
~~~

Every higher-order function also has a macro equivalent that accepts the function in-line rather than as a parameter:

~~~lisp
(funcall (bisect+ (x) (sin x)) 1 5)

3.1415926814079285d0
T
25

(funcall (bisect+ (x :tolerance 5d-12 :limit 100) (sin x)) 1 5)

3.1415926535919425d0
T
38
~~~

## Helper Functions

Since every closure returns 3 values (result, success/fail, and iteration count), some anaphoric macros
are included that can be used to make the multiple return values easier to work with:

~~~lisp
(if_ (funcall (bisect+ (x) (sin x)) 1 5)
     (format t "converged on ~A after ~A iterations.~%" _ @)
     (format t "unable to converge. Gave up after ~A iterations.~%" @))

converged on 3.1415926814079285d0 after 25 iterations.
NIL
~~~

The macros bind "_" to the numerical result and "@" to the number of iterations. The macros will behave like
their namesakes if used outside of this library:

~~~lisp
(if_ t 1 2)

1

(if_ nil 1 2)

2
~~~

## API Reference

Coming soon.

## License

MIT
