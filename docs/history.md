# Amplify History

AmplifyP is probably the 5th re-write of
[Bill Engel](https://genetics.wisc.edu/staff/engels-william/)'s the Amplify
software.

The earliest mentioning of Amplify that I could find is a
[Usenet announcement](https://groups.google.com/g/bionet.software/c/jK-a1gVUCV0).
The date was on 8 May 1992. Bill Engels probably uploaded it to
[Info-Mac](https://archive.org/details/info-mac-archive) earlier than that,
however I cannot find a reference. I also cannot locate the original Amplify
program (`amplify-10.hqx`). There are some evidences of
[early discussion](https://groups.google.com/g/bionet.software/c/ivIzidC7wpM) on
the original Amplify program. However they are hard to find, due to how long ago
it was. Amplify was published in a journal [^1] in November 1993.

[Amplify 3.2](https://github.com/wrengels/Amplify) was written in RealBasic. I
have had difficulty with finding the toolchain for compiling it for modern
operation system. The compiled binary can only be run on 32-bit Intel Mac. The
last version of macOS that could run this program was macOS
[Mojave](https://en.wikipedia.org/wiki/MacOS_Mojave), which was released on 24
Setepmber 2018. The end of life date for that version of macOS was in October
2021\.

[Amplify 4](https://github.com/wrengels/Amplify4) was written in a old dialect
of Swift. The fact that it was written in Swift means that Amplify 4 only works
on macOS. The fact that it was written in an *old dialect* of Swift means that
the source code cannot be re-compiled under modern Xcode. While I worked on
enabling Amplify 4 to
[support longer primer](https://github.com/wrengels/Amplify4/pull/1), I did
manage to recompile Amplify 4. However I do not recall the version of macOS and
Xcode that were needed.

I decided that it was probably for the best to rewrite Amplify 4 in Python,
because:

1. Python is a one of the most widely taught and used programming language.
   Choosing a popular language ensures the long-term survival of the project.
2. Python is cross-platform. This enables non-Mac user to use Amplify.

[^1]: Engels W(1993),
    [Contributing software to the internet: the Amplify program](<https://www.cell.com/trends/biochemical-sciences/abstract/0968-0004(93)90148-G>)
    *Trends in Biochemical Sciences*, 18, 448-450
