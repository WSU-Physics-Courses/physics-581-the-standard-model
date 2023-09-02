(sec:sylabus)=
# Syllabus

## Course Information

- **Instructor(s):** Michael McNeil Forbes [`m.forbes+581@wsu.edu`](mailto:m.forbes+581@wsu.edu)
- **Course Assistants:** 
- **Office:** Webster 947F
- **Office Hours:** TBD
- **Course Homepage:** https://schedules.wsu.edu/List/Pullman/20233/Phys/581/01
- **Class Number:** 581
- **Title:** Phys 581: The Standard Model of Particle Physics
- **Credits:** 3
- **Recommended Preparation**: 
- **Meeting Time and Location:** MW, 4:10pm - 4:00pm, [Webster 941](https://sah-archipedia.org/buildings/WA-01-075-0008-13), Washington State University (WSU), Pullman, WA
- **Grading:** Grade based on assignments and project presentation.

```{contents}
```

### Prerequisites
<!-- Only required for undergraduate courses, but is useful for graduate courses -->

This course is intended for a broad audience.  As such, the only formal background
assumed is a strong background in fundamental mathematical techniques, including complex
analysis, calculus, differential equations, and linear algebra as described in
{ref}`sec:linear-algebra`.

Familiarity with the core physics areas of classical mechanics, quantum mechanics,
statistical mechanics, and electromagnetism would be helpful, but we will try to
supplement any missing knowledge in the course.  (Speak up if you have a gap!)

### Textbooks and Resources

#### Required
:::{margin}
ISBN:9780691223483
:::
[Donoghue and Sorbo: "A Prelude to Quantum Field Theory" (2022)][Donoghue:2022]
: {cite:p}`Donoghue:2022` This is a principle textbook for the course.  It contains a
  fairly accessible introduction to the types of manipulations needed to do quantum
  field theory calculations.  Unfortunately, the eBook and required reader are pretty
  awful, so I highly recommend purchasing a paper copy.  Occasionally these go on sale
  for half-price.

:::{margin}
EISBN:9781315170626
:::
[Langacker: "The Standard Model and Beyond (2017)][Langacker:2017]
: {cite:p}`Langacker:2017` This is a principle textbook for the course.  It contains a
  complete and fairly straightforward introduction to the form of The Standard Model, as
  well as a nice discussion of possible extensions.  The CRC Press has made this [open
  access][Langacker:2017], making this an economical option.  Some of this will be a bit
  tough going if you don't have the background, but my hope is that you can follow the
  discussion and use it as a reference once you complete the course.

:::{margin}
ISBN:9780367091729 EISBN:9780429499210
:::
[Georgi: "Lie Algebras in Particle Physics" (2019)][Georgi:2019]
: {cite:p}`Georgi:2019` This is a fun introduction to group theory and Lie algebras with
  a focus on building the Standard Model.  The CRC Press has made this [open
  access][Georgi:2019].
  
['t Hooft: "The Conceptual Basis of Quantum Field Theory"][t-Hooft:2016]
: {cite:p}`t-Hooft:2016` This is a set of notes by Gerard 't Hooft that briefly
  discusses some of the conceptual foundations and issues with QFT and the Standard
  Model. 't Hooft has one of the most beautiful and deep understandings of field theory,
  and does not shy away from some of the difficult aspects of the field.  This will
  probably be hard going, but should give you a flavor for some of the subtle and deep
  aspects of the field.

[Donoghue:2022]: <https://press.princeton.edu/books/paperback/9780691223483/a-prelude-to-quantum-field-theory>
[Langacker:2017]: <https://www.taylorfrancis.com/books/oa-mono/10.1201/b22175/standard-model-beyond-paul-langacker>
[t-Hooft:2016]: <https://webspace.science.uu.nl/~hooft101/lectures/basisqft.pdf>
[Georgi:2019]: <https://doi.org/10.1201/9780429499210>

#### Additional Resources

Here are some additional optional resources.  Although they are not official texts for
the course, I will often refer to sections of them.  Feel free to ask questions about
these in class.

* {cite:p}`Zee:2010`: "Quantum Field Theory in a Nutshell".  This is a fun book to read
  by one of the masters of the field.  However, Zee presents things in such a way that
  it is easy for readers to feel they understand more than they do.  This book should be
  read in conjunction with a more careful text, and used for perspective.
* {cite:p}`Zee:2023`: "Quantum Field Theory, As Simply As Possible".  Bedtime reading
  with minimal equations.  Another good complement for a more serious text to help you
  gain perspective (if you enjoy Zee's style).
* {cite:p}`Donoghue:2014`: "Dynamics of the Standard} Model".  This is a serious
  and complete presentation of the Standard Model, similar to {cite:p}`Langacker:2017`,
  but more complete, and, consequently, more difficult to understand.  It is a great
  reference, but I fear will be very difficult for most students to understand.
* {cite:p}`Lancaster:2014`: "Quantum Field Theory for the Gifted Amateur".  A fairly
  easy textbook to start with if you want to learn QFT more thoroughly.
* {cite:p}`Maggiore:2005`: "A Modern Introduction to Quantum Field Theory".  Another
  good starting textbook that develops QFT technology starting from symmetries.
* {cite:p}`Nastase:2019`: "Classical Field Theory".  A comprehensive presentation of
  classical field theory – i.e. electromagnetism, gravitation, etc.  This is a good
  place to refresh your knowledge of relativity.

* {cite:p}`Cornwell:1997`: "Group Theory for Physics" A very fast and thorough
  introduction to group theory and its applications to physics.  Skips over proofs (with
  references) to emphasize the physics.
* {cite:p}`Boyd:1999` "The Devil's Invention: Asymptotic, Superasymptotic and
  Hyperasymptotic Series", and {cite:p}`Bender:1999` "Advanced Mathematical Methods for
  Scientists and Engineers I: Asymptotic Methods and Perturbation Theory": These provide
  technical details about asymptotic series of the type common in QFT. 
* {cite:p}`Huang:2013` A fun high-level discussion of renormalization.  Fig. 13 provides
  the basis for my general picture of physics.

Additional readings and references will be provided as needed.  Please see
{ref}`sec:readings` for details.
Details and further resources will also be included on the lecture pages on the
[Canvas](https://wsu.instructure.com/courses/1688766) server.



:::{margin}
In this course, the term "Standard Model" refers to the "Standard Model of Particle
Physics" in contradistinction to the "Standard Model of Cosmology" which includes the
components of inflation, the big bang, dark matter and dark energy (sometimes called
ΛCDM - for dark energy expressed as a cosmological constant "Λ" plus cold dark matter),
etc.

These two "Standard Models" are the most accurate and perplexing physical models we
have. They are perplexing in that we know they are wrong, but have no direct evidence
about how they breakdown.
:::
### Student Learning Outcomes

* **Dimensional Analysis**: All physical quantities are formed of four
  fundamental dimensions: **mass**, **time**, **distance**, and **charge**.  This
  affords a tremendously powerful technique called dimensional analysis that allows
  physicists to quickly estimate the form of an answer before doing any calculation.
  Dimensional analysis also allows one to check a calculation **or** to simplify a
  calculation by setting up to three independent quantities equal to unity.
  
  In discussions of the Standard Model, we often set $\hbar = c = 1$, thereby making
  equivalent energy, inverse time, and inverse length.

* **Symmetries and Conservation Laws**: Symmetry plays a deep role in the formulation of
  physical theories.  Through [Noether's theorem][], continuous symmetries lead to
  conservation laws.  For example:
  * Translation invariance ⟹ conservation of linear momentum;
  * Rotational invariance ⟹ conservation of angular momentum;
  * Time-translation invariance ⟹ conservation of energy.
  
  In the Standard Model, various other symmetries associated with Lie groups, give rise
  to various conserved **charges** that play a fundamental role.  (E.g. electric charge,
  baryon number, lepton number).  Symmetries also greatly restrict the allowed terms,
  leading to a complete specification of the model on a few simple properties.

[Noether's theorem]: <https://en.wikipedia.org/wiki/Noether's_theorem>

* Understand the ingredients and basic structure of the Standard Model, and how these
  come together to explain matter and its interactions.  This includes the following:
  * Quantum Field Theory (QFT) as the framework.
  * Fundamental space-time symmetries:
    * Relativistic Poincare invariance (Lorentz invariance + translations).
    * Non-relativistic limit of Galilean invariance (Rotations + boosts + translations).
    * CPT: Charge conjugation times parity times time-reversal invariance.  (Relates
      particles to anti-particles).
  * Representation of symmetries with groups.
  * Organization of matter in terms of spinors -- representations of the Poincare group.
    * Spin 0 fields (scalar bosons): the **Higgs** field.
    * Spin 1/2 fields (spinor fermions): **Electron**, **muon**, and **tauon**;
      Corresponding **Neutrinos**; **Quarks**.
    * Spin 1 fields (vector bosons): Mediators of force (interactions).  These vector
      bosons must be massless to avoid negative energies, and this ensured by imposing
      gauge symmetries  The following gauge bosons exist in the standard model:
      * $U(1)$: **Photon** (electromagnetic interactions)
      * $SU(2)$: **$W$ and $Z$ bosons** (weak interactions)
      * $SU(3)$: **Gluons** (strong interaction)

      *Note: The Higgs mechanism complicates the relationship betweem the bare gauge
      bosons that form the $SU(3)\times SU(2) \times U(1)$ gauge group of the standard
      model.*
    * Spin 2 fields: **Graviton**. Know that gravity fits within the framework of the
      Standard Model, but fails to be renormalizable (see below), rendering such a
      theory of gravity incomplete and lacks predictive power.
  * Principle of mass dimension and renormalizability.
* Understand how to combine these ingredients to form the Standard Model as the most
  general quantum field theory with these matter and gauge fields that satisfies all
  known symmetry constraints and includes only renormalizable terms.
* **QED**: Quantum electrodynamics consists of the QFT containing the electron and the
  photon.  Including atomic nuclei (made of protons and neutrons), this theory describes
  virtually all of the world we interact with daily, including most physical systems,
  chemistry, and biology.  This is the lowest energy scale of the Standard Model,
  dealing with atomic energies (eV to keV) and lengths scales on the order of
  angstroms (1Å=$10^{-10}$m) (think of the ground state of a Hydrogen atom etc.).
* **EW**: Adding to QED the $W$ and $Z$ bosons, as well as massive **mesons**
  (**pions**, **kaons**, and **tauons**) one obtains electroweak theory which includes
  phenomena such as radioactive decay: something we experience, but not on a daily basis
  (I hope!).  These occur at energies in the MeV range and lengths scales of picometers (1pm$=10^{-12}$m).
* **QCD**: Quantum chromodynamics includes the quarks and gluons as degrees of freedom.
  From these, we obtain the proton, neutron, and mesons included in the previous
  low-energy theories.  These interactions have energies of GeVs on the scale of a
  fermtometer (sometimes called a "fermi") 1fm$=10^{-15}$m.
* **Asymptotic Freedom**: QCD alone has another special feature -- the interactions
  become weaker and weaker as one goes to higher energy.  This means that as one probes
  with higher energies (equivalently, shorter length scales), the theory becomes
  non-interacting.  This introduces the tantalizing possibility that QCD is a correct
  theory "all the way down".  This is in contrast with QED, where the interactions get
  stronger at short distance, ultimately leading to a fundamental divergence (sometimes
  called the "Landau Pole") that tells us the theory cannot be complete mathematically.
* **Renormalization and Regularization:** From this description you see that the theory
  looks different at different energy scales.  This is the idea of renormalization


These ingredients and principles require many interesting and subtle mathematical
tools.  A good portion of this course will be devoted to explaining these tools, which
might have applications elsewhere.

* **Perturbation Theory and Feynman Diagrams:**
* **Path Integrals and Generating Functions:**
* **Green's functions:**
* **Asymptotic Series:**
* **Complex Analysis:**

### Expectations for Student Effort

For each hour of lecture equivalent, all students should expect to have a minimum of one
hour of work outside class.  All students are expected to keep up with the readings
assigned in class, asking questions through the Perusall/Hypothes.is forums, complete
homework on time, and prepare their projects/presentations.

### Assessment and Grading Policy

There are two options for obtaining a grade in this course.

* Complete all the assignments, obtaining a grade according to the following scale.
* Keep and submit a **class notebook** wherein you document your attempts to learn the
  material.  This notebook should contain at a minimum, dates and times when you work,
  state your learning goals, and then summarize what worked and what did not.  I expect
  that when you run into difficulties (as documented in your notebook), that you bring
  these issues to my attention in class, on the forums, or at office hours.  If you
  document a good attempt to learn the material covered in the assignments, then you
  will get an A for that assignment, even if you cannot successfully complete it.
  
  Your **class notebook** can be kept electronically on CoCalc, or you can submit scans
  of your physical notebook.
  
  :::{margin}
  Establishing precedence with an electronic notebook is a challenge.  One way to do
  this would be to publish regularly to a dated public repository such as [GitLab][].  A
  more secure approach might be to use and establish a
  [Blockchain](https://en.wikipedia.org/wiki/Blockchain).  This would be an interesting
  topic for discussion.
  :::
  You should treat this as a **lab notebook**, meaning that you should not edit it.  It
  should contain a chronological record of your progress through course, as an
  experimentalist would maintain chronological records of the experimental process.
  Such a notebook is extremely valuable.  For example, should you discover something, a
  notebook can be used to establish academic precedence.

The final grade will be converted to a letter grade using the following scale: 

| Percentage P       | Grade |
| ------------------ | ----- |
| 90.0% ≤ P          | A     |
| 85.0% ≤ P \< 90.0% | A-    |
| 80.0% ≤ P \< 85.0% | B+    |
| 75.0% ≤ P \< 80.0% | B     |
| 70.0% ≤ P \< 75.0% | B-    |
| 65.0% ≤ P \< 70.0% | C+    |
| 60.0% ≤ P \< 65.0% | C     |
| 55.0% ≤ P \< 60.0% | C-    |
| 50.0% ≤ P \< 55.0% | D+    |
| 40.0% ≤ P \< 50.0% | D     |
| P \< 40.0%         | F     |

::::{margin}
:::{warning}
Customize this section for your course!
:::
::::
### Attendance and Make-up Policy 

While there is no strict attendance policy, students are expected attend an participate
in classroom activities and discussion. Students who miss class are expected to cover
the missed material on their own, e.g. by borrowing their classmates notes, reviewing
recorded lectures (if available), etc.

### Course Timeline

As this is the first time we are offering this course as an iSciMath course, the
schedule well be extremely fluid in response to the experience of the class.

<!-- 16 Weeks -->

:::{margin}
*Weeks 1-3*
:::
**The Basis of Physics**
* Classical Mechanics
* Quantum Mechanics
* What is a Quantum Field Theory
* Zee's Baby Problem

:::{margin}
*Week 4*
:::
:::{margin}
*Week 5*
:::
:::{margin}
*Week 7*
:::
:::{margin}
*Week 8*
:::
:::{margin}
*Week 9*
:::
:::{margin}
*Week 10*
:::
:::{margin}
*Week 11*
:::
:::{margin}
*Week 12*
:::
:::{margin}
*Week 13*
:::
:::{margin}
*Week 14*
:::
*Thanksgiving Break -- No Classes* 

:::{margin}
*Weeks 15-16*
:::

**Course Review and Future Directions**

## Other Information

### Policy for the Use of Large Language Models (LLMs) or Generative AI in Physics Courses

The use of LLMs or Generative AI such as Chat-GPT is becoming prevalent, both in education and in industry.  As such, we believe that it is important for students to recognize the capabilities and inherent limitations of these tools, and use them appropriately.

To this end, **please submit 4 examples of your own devising:**
* Two of which demonstrate the phenomena of "hallucination" -- Attempt to use the tool to learn something you know to be true, and catch it making plausible sounding falsehoods.
* Two of which demonstrate something useful (often the end of a process of debugging and correcting the AI).

Note: one can find plenty of examples online of both cases.  Use these to better understand the capabilities and limitations of the AIs, but for your submission, please find your own example using things you know to be true. *If you are in multiple courses, you may submit the same four examples for each class, but are encouraged to tailor your examples to the course.*

Being able to independently establish the veracity of information returned by a search, an AI, or indeed any publication, is a critical skill for a scientist.  **If you are the type of employee who can use tools like ChatGPT to write prose, code etc., but not accurately validate the results, then you are exactly the type of employee that AI will be able to replace.** 

Any use of Generative AI or similar tools for submitted work must include:
1. **A complete description of the tool.** (E.g. *"ChatGPT Version 3.5 via CoCalc's interface"* or *Chat-GPT 4 through Bing AI using the Edge browser"*, etc.)
2. **A complete record of the queries issued and response provided.**  (This should be provided as an attachment, appendices, or supplement.)
3. **An attribution statement consistent with the following:**
   *“The author generated this <text/code/etc.> in part with <GPT-3, OpenAI’s large-scale language-generation model/etc.> as documented in appendix <1>. Upon generating the draft response, the author reviewed, edited, and revised the response to their own liking and takes ultimate responsibility for the content.”*

### Policy for the Use of Large Language Models (LLMs) or Generative AI in Physics Courses

The use of LLMs or Generative AI such as Chat-GPT is becoming prevalent, both in
education and in industry.  As such, we believe that it is important for students to
recognize the capabilities and inherent limitations of these tools, and use them
appropriately.

To this end, **please submit 4 examples of your own devising:**
* Two of which demonstrate the phenomena of "hallucination" -- Attempt to use the tool
  to learn something you know to be true, and catch it making plausible sounding
  falsehoods.
* Two of which demonstrate something useful (often the end of a process of debugging and
  correcting the AI).

Note: one can find plenty of examples online of both cases.  Use these to better
understand the capabilities and limitations of the AIs, but for your submission, please
find your own example using things you know to be true. *If you are in multiple courses,
you may submit the same four examples for each class, but are encouraged to tailor your
examples to the course.*

Being able to independently establish the veracity of information returned by a search,
an AI, or indeed any publication, is a critical skill for a scientist.  **If you are the
type of employee who can use tools like ChatGPT to write prose, code etc., but not
accurately validate the results, then you are exactly the type of employee that AI will
be able to replace.** 

Any use of Generative AI or similar tools for submitted work must include:
1. **A complete description of the tool.** (E.g. *"ChatGPT Version 3.5 via CoCalc's
   interface"* or *Chat-GPT 4 through Bing AI using the Edge browser"*, etc.)
2. **A complete record of the queries issued and response provided.**  (This should be
   provided as an attachment, appendices, or supplement.)
3. **An attribution statement consistent with the following:**
   *“The author generated this <text/code/etc.> in part with <GPT-3, OpenAI’s
   large-scale language-generation model/etc.> as documented in appendix <1>. Upon
   generating the draft response, the author reviewed, edited, and revised the response
   to their own liking and takes ultimate responsibility for the content.”*

<!-- ### COVID-19 Statement -->
<!-- Per the proclamation of Governor Inslee on August 18, 2021, **masks that cover both the -->
<!-- nose and mouth must be worn by all people over the age of five while indoors in public -->
<!-- spaces.**  This includes all WSU owned and operated facilities. The state-wide mask mandate -->
<!-- goes into effect on Monday, August 23, 2021, and will be effective until further -->
<!-- notice.  -->
 
<!-- Public health directives may be adjusted throughout the year to respond to the evolving -->
<!-- COVID-19 pandemic. Directives may include, but are not limited to, compliance with WSU’s -->
<!-- COVID-19 vaccination policy, wearing a cloth face covering, physically distancing, and -->
<!-- sanitizing common-use spaces.  All current COVID-19 related university policies and -->
<!-- public health directives are located at -->
<!-- [https://wsu.edu/covid-19/](https://wsu.edu/covid-19/).  Students who choose not to -->
<!-- comply with these directives may be required to leave the classroom; in egregious or -->
<!-- repetitive cases, student non-compliance may be referred to the Center for Community -->
<!-- Standards for action under the Standards of Conduct for Students. -->

### Academic Integrity

You are responsible for reading WSU’s [Academic Integrity Policy][], which is based on
[Washington State law][]. If you cheat in your work in this class you will: 

* Fail the course.
* Be reported to the [Center for Community Standards][].
* Have the right to appeal the instructor's decision.
* Not be able to drop the course or withdraw from the course until the appeals process
  is finished.

If you have any questions about what you can and cannot do in this course, ask your instructor.

If you want to ask for a change in the instructor's decision about academic integrity,
use the form at the [Center for Community Standards][] website. You must submit this
request within 21 calendar days of the decision.

[Academic Integrity Policy]: 
  <https://communitystandards.wsu.edu/policies-and-reporting/academic-integrity-policy/>
[Washington State law]: <https://apps.leg.wa.gov/wac/default.aspx?cite=504-26-202>
[Center for Community Standards]: <https://communitystandards.wsu.edu/>

<!-- Academic integrity is the cornerstone of higher education.  As such, all members of the -->
<!-- university community share responsibility for maintaining and promoting the principles -->
<!-- of integrity in all activities, including academic integrity and honest -->
<!-- scholarship. Academic integrity will be strongly enforced in this course.  Students who -->
<!-- violate WSU's Academic Integrity Policy (identified in Washington Administrative Code -->
<!-- (WAC) [WAC 504-26-010(4)][wac 504-26-010(4)] and -404) will fail the course, will not -->
<!-- have the option to withdraw from the course pending an appeal, and will be Graduate: -->
<!-- 6300, 26300 reported to the Office of Student Conduct. -->

<!-- Cheating includes, but is not limited to, plagiarism and unauthorized collaboration as -->
<!-- defined in the Standards of Conduct for Students, [WAC 504-26-010(4)][wac -->
<!-- 504-26-010(4)]. You need to read and understand all of the [definitions of -->
<!-- cheating][definitions of cheating].  If you have any questions about what is and is not -->
<!-- allowed in this course, you should ask course instructors before proceeding. -->

<!-- If you wish to appeal a faculty member's decision relating to academic integrity, please -->
<!-- use the form available at \_communitystandards.wsu.edu. Make sure you submit your appeal -->
<!-- within 21 calendar days of the faculty member's decision. -->

<!-- Academic dishonesty, including all forms of cheating, plagiarism, and fabrication, is -->
<!-- prohibited. Violations of the academic standards for the lecture or lab, or the -->
<!-- Washington Administrative Code on academic integrity -->

### University Syllabus
<!-- Required as of Fall 2023 -->

Students are responsible for reading and understanding all university-wide policies and
resources pertaining to all courses (for instance: accommodations, care resources,
policies on discrimination or harassment), which can be found in the [university
syllabus][].

[university syllabus]: <https://syllabus.wsu.edu/university-syllabus/>

### Students with Disabilities

Reasonable accommodations are available for students with a documented
disability. If you have a disability and need accommodations to fully
participate in this class, please either visit or call the Access
Center at (Washington Building 217, Phone: 509-335-3417, E-mail:
<mailto:Access.Center@wsu.edu>, URL: <https://accesscenter.wsu.edu>) to schedule
an appointment with an Access Advisor. All accommodations MUST be
approved through the Access Center. For more information contact a
Disability Specialist on your home campus.

### Campus Safety

Classroom and campus safety are of paramount importance at Washington
State University, and are the shared responsibility of the entire
campus population. WSU urges students to follow the “Alert, Assess,
Act,” protocol for all types of emergencies and the “[Run, Hide, Fight]”
response for an active shooter incident. Remain ALERT (through direct
observation or emergency notification), ASSESS your specific
situation, and ACT in the most appropriate way to assure your own
safety (and the safety of others if you are able).

Please sign up for emergency alerts on your account at MyWSU. For more
information on this subject, campus safety, and related topics, please
view the FBI’s [Run, Hide, Fight] video and visit [the WSU safety
portal][the wsu safety portal].

### Students in Crisis - Pullman Resources 

If you or someone you know is in immediate danger, DIAL 911 FIRST! 

* Student Care Network: https://studentcare.wsu.edu/
* Cougar Transit: 978 267-7233 
* WSU Counseling and Psychological Services (CAPS): 509 335-2159 
* Suicide Prevention Hotline: 800 273-8255 
* Crisis Text Line: Text HOME to 741741 
* WSU Police: 509 335-8548 
* Pullman Police (Non-Emergency): 509 332-2521 
* WSU Office of Civil Rights Compliance & Investigation: 509 335-8288 
* Alternatives to Violence on the Palouse: 877 334-2887 
* Pullman 24-Hour Crisis Line: 509 334-1133 

[communitystandards.wsu.edu]: https://communitystandards.wsu.edu/
[definitions of cheating]: https://apps.leg.wa.gov/WAC/default.aspx?cite=504-26-010
[run, hide, fight]: https://oem.wsu.edu/emergency-procedures/active-shooter/
[the wsu safety portal]: https://oem.wsu.edu/about-us/
[wac 504-26-010(4)]: https://apps.leg.wa.gov/WAC/default.aspx?cite=504-26
[SSH]: <https://en.wikipedia.org/wiki/Secure_Shell> "SSH on Wikipedia"
[CoCalc]: <https://cocalc.com> "CoCalc: Collaborative Calculation and Data Science"
[GitLab]: <https://gitlab.com> "GitLab"
[GitHub]: <https://github.com> "GitHub"
[Git]: <https://git-scm.com> "Git"
[Anki]: <https://apps.ankiweb.net/> "Powerful, intelligent flash cards."


[Official Course Repository]: <https://gitlab.com/wsu-courses/physics-581-the-standard-model> "Official Course Repository hosted on GitLab"
[Shared CoCalc Project]: <https://cocalc.com/4578ccc8-55cf-413d-9e4c-82c69e30121e/> "Shared CoCalc Project"
[WSU Courses CoCalc project]: <https://cocalc.com/projects/c31d20a3-b0af-4bf7-a951-aa93a64395f6>


