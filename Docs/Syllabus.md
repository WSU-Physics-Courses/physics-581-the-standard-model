(sec:sylabus)=
# Syllabus

::::{margin}
:::{warning}
Customize this section for your course!
:::
::::
::::{admonition} Instructor Notes
This syllabus contains sample content from the WSU course **Physics 455/555: Quantum
Technologies and Computation** as offered Fall 2022.  These sections need to be replaced
with appropriate material for your course.  Material that needs to be customized
has a "Customize" margin note.  In the code look for: 

    ::::{margin}
    :::{warning}
    Customize this section for your course!
    :::
    ::::

Once the material is suitably customize, just remove the margin note.
::::

::::{margin}
:::{warning}
Customize this section for your course!  Most of this should come from the cookiecutter
template `.cookiecutter.yaml`, but it might need some attention.

:::
::::
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

:::{admonition} Customize for Course
This course is intended for a broad audience.  As such, the only formal background
assumed is a strong background in linear algebra as described in {ref}`sec:linear-algebra`.

Familiarity with quantum mechanics, classical information theory, and the foundations of
computer science may be useful, but are not required.  The course will focus on
finite-dimensional systems where a deep understanding of complex vector spaces and
matrices will suffice.
:::

::::{margin}
:::{warning}
Customize this section for your course!
:::
::::
### Textbooks and Resources

#### Required
:::{margin}
ISBN:9781107002173, EISBN:9780511990809
:::
[Nielsen and Chuang: "Quantum Computation and Quantum Information" (2010)][Nielsen:2010 (ProQuest)]
:  This is the principle textbook for the course by two of the founders of the field.
   It provides a reasonably accessible, but very thorough review of the field, replete
   with references and history.  It may be a bit dense on first reading, but lays out a
   complete self-contained foundation.  It is available through the [WSU
   Library][Nielsen:2010 (WSU)] electronically through [ProQuest][Nielsen:2010
   (ProQuest)] (requires WSU sign-in).

:::{margin}
The basic idea is that key concepts are presented as simple questions that you are
repeatedly asked.  If you remember the concepts, then the frequency at which you are
asked is reduced.  This is similar to using flash-cards, but with an automated delivery
schedule optimized by cognition research for maximal retention.

For a discussion, see [Augmenting Long-term Memory](http://augmentingcognition.com/ltm.html)
:::
[Quantum Country][] by Andy Matuschak and Michael Nielsen
:  A unique presentation of some of the key concepts in quantum computing and quantum
   mechanics by one of the authors of the principle textbook.  This is required reading
   for the course.  It is not a traditional textbook, nor is it a complete presentation
   of the material, but presents key concepts using a new mnemonic medium to help you
   remember the core concepts.  We use this mnemonic medium in our notes and encourage
   you to explore further using tools like [Anki][] to make notes for yourself (but see
   [Augmenting Long-term Memory][] -- there is an art to making good notes).


#### Additional Resources

[Preskill: "Quantum Computation"][Preskill (Physics/CS 219)]
:  Online notes from John Preskill.  These provide more theoretical foundation that
   complements {cite:p}`Nielsen:2010`.

Additional readings and references will be provided as needed.  Please see
{ref}`sec:readings` for details.
Details and further resources will also be included on the lecture pages on the
[Canvas](https://wsu.instructure.com/courses/1688766) server.


[Augmenting Long-term Memory]: <http://augmentingcognition.com/ltm.html> 
  "Michael Nielsen's discussion about how to use Anki effectively"
[Quantum Country]: <https://quantum.country/>
[Nielsen:2010 (ProQuest)]:
  <https://ntserver1.wsulibs.wsu.edu:2171/lib/wsu/detail.action?docID=647366>
[Nielsen:2010 (WSU)]:
  <https://searchit.libraries.wsu.edu/permalink/f/1jnr272/TN_cdi_askewsholts_vlebooks_9780511985249>
[Preskill (Physics/CS 219)]:
  <http://www.theory.caltech.edu/~preskill/ph229/>


::::{margin}
:::{warning}
Customize this section for your course!
:::
::::
### Student Learning Outcomes

**Physics 455/555:** By the end of this course, all students will:

1. Know the postulates of quantum mechanics and their consequences. 
2. Be able to analyze quantum circuits. 
3. Be familiar with the fundamental quantum algorithms. 
4. Be familiar with current technologies being explored for realizing quantum computing. 
5. Be able to use a quantum simulation or computing platform to implement quantum algorithms.
6. Be aware of the potential advantages offered by quantum technologies, but cognizant
   of the physical challenges and limitations.

**Physics 555:** Graduate students will additionally:

7. Be able to review the literature about a specific topic relevant to the course that
   is under active investigate, and present a critical summary of this topic to the
   class considering the previous outcomes.

::::{margin}
:::{warning}
Customize this section for your course!
:::
::::
### Expectations for Student Effort

For each hour of lecture equivalent, all students should expect to have a minimum of two
hours of work outside class.  All students are expected to keep up with the readings
assigned in class, asking questions through the Perusall/Hypothes.is forums, complete
homework on time, and prepare their projects/presentations.

::::{margin}
:::{warning}
Customize this section for your course!
:::
::::
### Assessment and Grading Policy

Assessment and Grading Policy 

**Physics 455/555:** Students will be assessed with 5 weekly assignments to ensure that
the content-based learning outcomes 1 through 4 are realized.  The first four
assignments will be due at the start of week following that when the material is
presented. The fifth assignment will have students complete exercises using a platform
like Qiskit (learning outcome 5).  These assignments will be worth 70% of the final
grade (15% each for assignments 1-4, 10% for assignment 5). Late assignments will not be
accepted unless prior arrangements are made with the instructor. The course will have no
examinations.

**Physics 455:** Undergraduates will be required to complete a quantum computation
project, worth the remaining 30% of the student’s grade, using one of the available
quantum simulation or computing platforms, and to report on this, discussing any
limitations imposed by real hardware, errors, decoherence, etc. (Learning outcomes 5 and
6). Projects are due at the start of the final week of class. Late submissions will not
be accepted unless prior arrangements are made with the instructor.

**Physics 555:** Graduate students must make a presentation, worth the remaining 30% of
the student’s grades, about an active research topic relevant to the course, subject to
instructor approval, and must demonstrate a review of the appropriate literature
relevant to this topic. (Learning outcome 7).  These presentations will be made in weeks
5-16 of the course, interspersed with related discussions and background material.
Presentations might take the form of a lecture or in-class presentation, or could
consist of a set of notes in the style of the course documentation.  Both must include
references.

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

::::{margin}
:::{warning}
Customize this section for your course!
:::
::::
### Course Timeline
<!-- 16 Weeks -->

:::{margin}
*Week 1 / Assignment 1*
:::
Linear Algebra
: - Eigenvectors and Eigenvalues
  - Hermitian

:::{margin}
*Week 2 / Assignment 2*
:::
Quantum Mechanics:
:  -  Postulates.
   -  Spin ½, Pauli matrices, Rotations.
   -  EPR pairs.
   -  Measurement
      -  POV Measurements
      -  von Neumann measurements (pointer states)
   - Entanglement Measures
   - Bell's Inequalities

:::{margin}
*Week 3 / Assignment 3*
:::
Quantum Computing Theory
:  -  No Cloning
   -  Universal gates.
   -  Communication
      -  Alice, Bob, and Eve
      -  Superdense coding
      -  Teleportation
   
:::{margin}
*Week 4 / Assignment 4*
:::

Quantum Circuits and Algorithms
:  - Grover, Deutsch, Schor
   - Adiabatic quantum computing
   - Quantum annealing


:::{margin}
*Weeks 5-13 / Assignment 5 (Qiskit) / 555 Presentations start*

*Week 7 / 455 Project Approval*

:::
Additional topics and presentations
: The remaining classes will contain a mix of lectures, guest lectures, and graduate
  student presentations about topics of interest to the class.  Exact content will be
  tailored to the interests of the current class and may include the following (which
  provides a partial list of potential topics for graduate student presentations):
  - Classical programming models and complexity
  - Quantum programming:  Basic gates, and programming models (universal computing)
  - Quantum Algorithms
  - Ion Trapping
  - Quantum Error Correction: Ancillary qubits
  - Quantum Complexity: What can be done with quantum computers that cannot be done with
    classical computers (presumably and provably)
  - Communication: Encryption protocols
  - Quantum Optics: Experimental quantum communication
  - Neutral atoms in lattices and optical tweezers
  - Atom interferometry 
  - Entanglement purification/Producing ground states
  - Superconducting qubits
  - NV centers
  - MRI quantum computing
  - Quantum Simulation

:::{margin}
*Week 14*
:::
*Thanksgiving Break -- No Classes* 

:::{margin}
*Weeks 15-16*
:::
Quantum Computing and Presentations
: **Physics 455:** Undergraduate students will run their computing projects on a quantum
  simulation or computing platform.
  
  **Physics 555:** Remaining students will deliver their presentations to the class.

### Physics 555 Presentation Topics

-  Presentations should be 30-40 minutes to leave time for discussion.
-  Please create a folder for your topic in "Notes and Discussions", and upload any
   slides, notes, or other relevant material for your topic.
-  Please create a Notes document in your folder with the following:
   -  **Introduction:** Please provide a summary of the topic.
   -  **References:** Include any helpful references.
   -  **Questions:** Please curate questions and discussions in this section. Use this
      section to ask any questions you have while working on your material, and then to
      help answer classmate questions.  This should be a dynamic section where all
      members of the course participate.

## Other Information

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
