---
layout: workshop      # DON'T CHANGE THIS.
carpentry: "dc"    # what kind of Carpentry (must be either "lc" or "dc" or "swc")
venue: "RNA-Seq workshop - University of Amsterdam"        # brief name of host site without address (e.g., "Euphoric State University")
address: "Room G3.10, Science Park 904, 1098XH Amsterdam, The Netherlands"      
country: "nl"      # lowercase two-letter ISO country code such as "fr" (see https://en.wikipedia.org/wiki/ISO_3166-1)
language: "en"     # lowercase two-letter ISO language code such as "fr" (see https://en.wikipedia.org/wiki/ISO_639-1)
latlng: "52.355354, 4.956575"       # decimal latitude and longitude of workshop venue (e.g., "41.7901128,-87.6007318" - use https://www.latlong.net/)
humandate: "March 7, 2019"    # human-readable dates for the workshop (e.g., "Feb 17-18, 2020")
humantime: "13:00 - 17:00"    # human-readable times for the workshop (e.g., "9:00 am - 4:30 pm")
startdate: 2019-03-07      # machine-readable start date for the workshop in YYYY-MM-DD format like 2015-01-01
enddate: 2019-03-07        # machine-readable end date for the workshop in YYYY-MM-DD format like 2015-01-02
instructor: ["Ernest Aliche", "Tijs Bliek", "Marc Galland"] # boxed, comma-separated list of instructors' names as strings, like ["Kay McNulty", "Betty Jennings", "Betty Snyder"]
helper: ["Max van Hooren","Fred White"]     # boxed, comma-separated list of helpers' names, like ["Marlyn Wescoff", "Fran Bilas", "Ruth Lichterman"]
email: ["m.galland@uva.nl","bliek@uva.nl"]    # boxed, comma-separated list of contact email addresses for the host, lead instructor, or whoever else is handling questions, like ["marlyn.wescoff@example.org", "fran.bilas@example.org", "ruth.lichterman@example.org"]
collaborative_notes:      # optional: URL for the workshop collaborative notes, e.g. an Etherpad or Google Docs document
eventbrite: 56550248315          # optional: alphanumeric key for Eventbrite registration, e.g., "1234567890AB" (if Eventbrite is being used)
---

{% comment %} See instructions in the comments below for how to edit specific sections of this workshop template. {% endcomment %}

{% comment %}
  HEADER
  Edit the values in the block above to be appropriate for your workshop.
  If the value is not 'true', 'false', 'null', or a number, please use
  double quotation marks around the value, unless specified otherwise.
  And run 'make workshop-check' *before* committing to make sure that changes are good.
{% endcomment %}

{% comment %}
  EVENTBRITE
  This block includes the Eventbrite registration widget if
  'eventbrite' has been set in the header.  You can delete it if you
  are not using Eventbrite, or leave it in, since it will not be
  displayed if the 'eventbrite' field in the header is not set.
{% endcomment %}
{% if page.eventbrite %}
<iframe
  src="https://www.eventbrite.com/tickets-external?eid={{page.eventbrite}}&ref=etckt"
  frameborder="0"
  width="100%"
  height="280px"
  scrolling="auto">
</iframe>
{% endif %}

{% comment %}
<h4>This is the workshop template. Delete these lines and use it to customize your own website.
If you are running a self-organized workshop or have not put in a workshop request yet, please also fill in
<a href="{{site.amy_site}}/submit">this workshop request form</a> to let us know about your workshop
and our administrator may contact you if we need any extra information.</h4>
{% endcomment %}

<h2 id="general">General Information</h2>

{% comment %}
  ############### INTRODUCTION ###################
  Edit the general explanatory paragraph below if you want to change the pitch.
{% endcomment %}
{% if page.carpentry == "swc" %}
  {% include sc/intro.html %}
{% elsif page.carpentry == "dc" %}
  {% include dc/intro.html %}
{% elsif page.carpentry == "lc" %}
  {% include lc/intro.html %}
{% endif %}

{% comment %}
  ############# AUDIENCE ############
Explain who your audience is.  (In particular, tell readers if the workshop is only open to people from a particular institution.
{% endcomment %}

{% if page.carpentry == "swc" %}
  {% include sc/who.html %}
{% elsif page.carpentry == "dc" %}
  {% include dc/who.html %}
{% elsif page.carpentry == "lc" %}
  {% include lc/who.html %}
{% endif %}


{% comment %}
  ############### LOCATION #############
  This block displays the address and links to maps showing directions
  if the latitude and longitude of the workshop have been set.  You
  can use https://itouchmap.com/latlong.html to find the lat/long of an
  address.
{% endcomment %}

{% if page.latlng %}
<p id="where">
<strong>Where:</strong> The workshop will be located in room G3.10 ({{page.address}}). Find the exact locations on <a href="//www.openstreetmap.org/?mlat={{page.latlng | replace:',','&mlon='}}&zoom=16">OpenStreetMap</a> or <a href="//maps.google.com/maps?q={{page.latlng}}">Google Maps</a>.
</p>
{% endif %}


{% comment %}
  ############# REQUIREMENTS ##########
{% endcomment %}

<p id="requirements">
  <strong>Requirements:</strong> Participants must bring a laptop with a
  Linux, Mac or Windows operating system (not a tablet, Chromebook, etc.) that they have administrative privileges on. They should have a few specific software packages installed (listed <a href="#setup">below</a>). They are also required to abide by
  {% if page.carpentry == "swc" %}
  Software Carpentry's
  {% elsif page.carpentry == "dc" %}
  the following
  {% elsif page.carpentry == "lc" %}
  Library Carpentry's
  {% endif %}
  <a href="{{site.swc_site}}/conduct.html">Code of Conduct</a>.
</p>

<p><strong>Registration:</strong> All free for Green Life Sciences students and staff from the University of Amsterdam. If not all seats are claimed by the <a href="http://gls.uva.nl/">Green Life Sciences cluster</a>, then the remaining seats will be open for other UvA students and staff. Please get your ticket using the Eventbrite widget (see above).</p>

<p><strong>Organisers:</strong> This workshop is organized by the <a href="https://www.scienceparkstudygroup.info/">Science Park Study Group</a> members as part of its recurring activities. Please visit our website for other events that take place in Amsterdam (usually Science Park). If there are any questions get in touch and send us an email (see <a href="#contact">contact</a>) </p>

{% comment %}
##########  ACCESSIBILITY ###############
{% endcomment %}

<p id="accessibility">
  <strong>Accessibility:</strong> We are committed to making this workshop
  accessible to everybody.
  The workshop organizers have checked that:
</p>
<ul>
  <li>The room is wheelchair / scooter accessible.</li>
  <li>Accessible restrooms are available.</li>
</ul>

{% comment %}
<p>
  Materials will be provided in advance of the workshop and large-print handouts are available if needed by notifying the
  organizers in advance.  If we can help making learning easier for you (e.g. sign-language interpreters, lactation facilities) please get in touch (using contact details below) and we will attempt to provide them.
</p>
{% endcomment %}

{% comment %}
  CONTACT EMAIL ADDRESS

  Display the contact email address set in the configuration file.
{% endcomment %}
<p id="contact">
  <strong>Contact</strong>:
  Please email
  {% if page.email %}
    {% for email in page.email %}
      {% if forloop.last and page.email.size > 1 %}
        or
      {% else %}
        {% unless forloop.first %}
        ,
        {% endunless %}
      {% endif %}
      <a href='mailto:{{email}}'>{{email}}</a>
    {% endfor %}
  {% else %}
    to-be-announced
  {% endif %}
  for more information.
</p>

<hr/>

{% comment %}
  ######### SCHEDULE ###############
  Show the workshop's schedule.  Edit the items and times in the table
  to match your plans.  You may also want to change 'Day 1' and 'Day
  2' to be actual dates or days of the week.
{% endcomment %}
<h2 class="schedule-title">Afternoon schedule</h2>

{% if page.carpentry == "swc" %}
  {% include  python/schedule.html %}
{% elsif page.carpentry == "dc" %}
  {% include dc/schedule.html %}
{% elsif page.carpentry == "lc" %}
  {% include lc/schedule.html %}
{% endif %}


<hr/>

{% comment %}
  ########### SYLLABUS ###############
  Show what topics will be covered.

  1. If your workshop is R rather than Python, remove the comment
     around that section and put a comment around the Python section.
  2. Some workshops will delete SQL.
  3. Please make sure the list of topics is synchronized with what you
     intend to teach.
  4. You may need to move the div's with class="col-md-6" around inside
     the div's with class="row" to balance the multi-column layout.

  This is one of the places where people frequently make mistakes, so
  please preview your site before committing, and make sure to run
  'tools/check' as well.
{% endcomment %}
<h2 id="syllabus">Syllabus</h2>

{% if page.carpentry == "swc" %}
  {% include python/syllabus.html %}
{% elsif page.carpentry == "dc" %}
  {% include dc/syllabus.html %}
{% elsif page.carpentry == "lc" %}
  {% include lc/syllabus.html %}
{% endif %}

<br>
<hr/>

{% comment %}
############## SETUP #############
  SETUP
  Delete irrelevant sections from the setup instructions.  Each
  section is inside a 'div' without any classes to make the beginning
  and end easier to find.
  This is the other place where people frequently make mistakes, so
  please preview your site before committing, and make sure to run
  'tools/check' as well.
{% endcomment %}

<h2 id="setup">Setup</h2>

<h3> Required Files </h3>
You can download all files (fastq sequencing files, genome sequence in fasta format, genome annotation, etc.) <a href="">here</a>. Unzip/extract the files to create a <em>rnaseq</em> folder.


<h3> Softwares</h3>
<p>
  To participate to this
  {% if page.carpentry == "swc" %}
  Software Carpentry
  {% elsif page.carpentry == "dc" %}
  RNA-Seq workshop
  {% elsif page.carpentry == "lc" %}
  Library Carpentry
  {% endif %}
  workshop,
  you will need access to the softwares described below.
  In addition, you will need an up-to-date web browser.
</p>
<p>
  We maintain a list of common issues that occur during installation as a reference for instructors
  that may be useful on the
  <a href = "{{site.swc_github}}/workshop-template/wiki/Configuration-Problems-and-Solutions">Configuration Problems and Solutions wiki page</a>.
</p>

<div id="shell"> {% comment %} Start of 'shell' section. {% endcomment %}
  <h3>The Bash Shell</h3>

  <p>
    Bash is a commonly-used shell that gives you the power to do simple tasks more quickly.
  </p>

  <div class="row">
    <div class="col-md-4">
      <h4 id="shell-windows">Windows</h4>
      <a href="https://www.youtube.com/watch?v=339AEqk9c-8">Video Tutorial</a>
      <ol>
        <li>Download the Git for Windows <a href="https://git-for-windows.github.io/">installer</a>.</li>
        <li>Run the installer and follow the steps bellow:
          <ol>
            {% comment %} Git 2.18.0 Setup {% endcomment %}
            <li>
                Click on "Next" four times (two times if you've previously
                installed Git).  You don't need to change anything
                in the Information, location, components, and start menu screens.
            </li>
            <li>
                <strong>
                Select “Use the nano editor by default” and click on “Next”.
                </strong>
            </li>
            {% comment %} Adjusting your PATH environment {% endcomment %}
            <li>
                Keep "Use Git from the Windows Command Prompt" selected and click on "Next".
                If you forgot to do this programs that you need for the workshop will not work properly.
                If this happens rerun the installer and select the appropriate option.
            </li>
            {% comment %} Choosing the SSH executable {% endcomment %}
            <li>Click on "Next".</li>
            {% comment %} Configuring the line ending conversions {% endcomment %}
            <li>
                Keep "Checkout Windows-style, commit Unix-style line endings" selected and click on "Next".
            </li>
            {% comment %} Configuring the terminal emulator to use with Git Bash {% endcomment %}
            <li>
              <strong>
                Select "Use Windows' default console window" and click on "Next".
              </strong>
            </li>
            {% comment %} Configuring experimental performance tweaks {% endcomment %}
            <li>Click on "Install".</li>
            {% comment %} Installing {% endcomment %}
            {% comment %} Completing the Git Setup Wizard {% endcomment %}
            <li>Click on "Finish".</li>
          </ol>
        </li>
        <li>
          If your "HOME" environment variable is not set (or you don't know what this is):
          <ol>
            <li>Open command prompt (Open Start Menu then type <code>cmd</code> and press [Enter])</li>
            <li>
              Type the following line into the command prompt window exactly as shown:
              <p><code>setx HOME "%USERPROFILE%"</code></p>
            </li>
            <li>Press [Enter], you should see <code>SUCCESS: Specified value was saved.</code></li>
            <li>Quit command prompt by typing <code>exit</code> then pressing [Enter]</li>
          </ol>
        </li>
      </ol>
      <p>This will provide you with both Git and Bash in the Git Bash program.</p>
    </div>
    <div class="col-md-4">
      <h4 id="shell-macosx">macOS</h4>
      <p>
        The default shell in all versions of macOS is Bash, so no
        need to install anything.  You access Bash from the Terminal
        (found in
        <code>/Applications/Utilities</code>).
        See the Git installation <a href="https://www.youtube.com/watch?v=9LQhwETCdwY ">video tutorial</a>
        for an example on how to open the Terminal.
        You may want to keep
        Terminal in your dock for this workshop.
      </p>
    </div>
    <div class="col-md-4">
      <h4 id="shell-linux">Linux</h4>
      <p>
        The default shell is usually Bash, but if your
        machine is set up differently you can run it by opening a
        terminal and typing <code>bash</code>.  There is no need to
        install anything.
      </p>
    </div>
  </div>
</div> {% comment %} End of 'shell' section. {% endcomment %}



<div id="editor"> {% comment %} Start of 'editor' section. {% endcomment %}
  <h3>Text Editor</h3>

  <p>
    When you're writing code, it's nice to have a text editor that is
    optimized for writing code, with features like automatic
    color-coding of key words.  The default text editor on macOS and
    Linux is usually set to Vim, which is not famous for being
    intuitive.  If you accidentally find yourself stuck in it, try
    typing the escape key, followed by <code>:q!</code> (colon, lower-case 'q',
    exclamation mark), then hitting Return to return to the shell.
  </p>

  <div class="row">
    <div class="col-md-4">
      <h4 id="editor-windows">Windows</h4>
      <p>
        nano is a basic editor and the default that instructors use in the workshop.
        It is installed along with Git.
      </p>
      <p>
        Others editors that you can use are
        <a href="https://notepad-plus-plus.org/">Notepad++</a> or
        <a href="https://www.sublimetext.com/">Sublime Text</a>.
        <strong>Be aware that you must
          add its installation directory to your system path.</strong>
        Please ask your instructor to help you do this.
      </p>
    </div>
    <div class="col-md-4">
      <h4 id="editor-macosx">macOS</h4>
      <p>
        nano is a basic editor and the default that instructors use in the workshop.
        See the Git installation <a href="https://www.youtube.com/watch?v=9LQhwETCdwY ">video tutorial</a>
        for an example on how to open nano.
        It should be pre-installed.
      </p>
      <p>
        Others editors that you can use are
        <a href="https://www.barebones.com/products/textwrangler/">Text Wrangler</a> or
        <a href="https://www.sublimetext.com/">Sublime Text</a>.
      </p>
    </div>
    <div class="col-md-4">
      <h4 id="editor-linux">Linux</h4>
      <p>
        nano is a basic editor and the default that instructors use in the workshop.
        It should be pre-installed.
      </p>
      <p>
        Others editors that you can use are
        <a href="https://wiki.gnome.org/Apps/Gedit">Gedit</a>,
        <a href="https://kate-editor.org/">Kate</a> or
        <a href="https://www.sublimetext.com/">Sublime Text</a>.
      </p>
    </div>
  </div>
</div> {% comment %} End of 'editor' section. {% endcomment %}

<h3> Credits</h3>
<p>
This material was built with heavy inspiration from the <a href="https://datacarpentry.org/genomics-workshop/">Data Carpentry Genomics workshop</a> developped under the supervision of the <a href="https://carpentries.org/">Carpentries Foundation</a>.  
We also borrowed material from a course developed by the <a href="https://hbctraining.github.io/Intro-to-rnaseq-hpc-O2/schedule/">Harvard Chan Bioinformatics Core</a>.
</p>
