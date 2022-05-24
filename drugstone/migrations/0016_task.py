# Generated by Django 3.0.5 on 2020-04-02 16:54

from django.db import migrations, models


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0015_auto_20200401_1717'),
    ]

    operations = [
        migrations.CreateModel(
            name='Task',
            fields=[
                ('id', models.AutoField(auto_created=True, primary_key=True, serialize=False, verbose_name='ID')),
                ('token', models.CharField(max_length=32, unique=True)),
                ('algorithm', models.CharField(max_length=128)),
                ('parameters', models.TextField()),
                ('progress', models.FloatField(default=0.0)),
                ('started_at', models.DateTimeField(null=True)),
                ('finished_at', models.DateTimeField(null=True)),
                ('worker_id', models.IntegerField(null=True)),
                ('done', models.BooleanField(default=False)),
                ('result', models.TextField(null=True)),
            ],
        ),
    ]